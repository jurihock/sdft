"""
Copyright (c) 2025 Juergen Hock

SPDX-License-Identifier: MIT

Sliding DFT implementation according to [1] combined with [2].

[1] Rick Lyons
    An Efficient Full-Band Sliding DFT Spectrum Analyzer
    DSPRelated.com (2011)
    https://www.dsprelated.com/showarticle/1396.php

[2] Russell Bradford and Richard Dobson and John ffitch
    Sliding is Smoother than Jumping
    International Computer Music Conference (2005)
    http://hdl.handle.net/2027/spo.bbp2372.2005.086

Source: https://github.com/jurihock/sdft
"""


import numba
import numpy


class SDFT:
    """
    Sliding Discrete Fourier Transform (SDFT).
    """

    def __init__(self, dftsize, window='hann', latency=1):
        """
        Create a new SDFT plan.

        Parameters
        ----------
        dftsize : int
            Desired number of DFT bins.
        window : str, optional
            Analysis window type (boxcar, hann, hamming or blackman).
        latency : float, optional
            Synthesis latency factor between 0 and 1.
            The default value 1 corresponds to the highest latency and best possible SNR.
            A smaller value decreases both latency and SNR, but also increases the workload.
        """

        if dftsize % 2:                          # odd dftsize
            fullsize = (dftsize - 1) * 2         # even fullsize
            assert dftsize == (fullsize / 2) + 1
        else:                                    # even dftsize
            fullsize = dftsize * 2 - 1           # odd fullsize
            assert dftsize == (fullsize + 1) / 2

        self.odd = fullsize % 2
        self.even = not self.odd

        self.size = dftsize
        self.window = window
        self.latency = latency

        self.buffer = numpy.zeros(dftsize, complex)

        self.twiddles_analysis = numpy.exp(+2j * numpy.pi * numpy.arange(dftsize) / fullsize)
        self.twiddles_synthesis = numpy.exp(-1j * numpy.pi * numpy.arange(dftsize) * latency)

        if self.latency == 1:

            # circular shift in time domain or multiplication of each dft bin by (-1)**n
            self.twiddles_synthesis = numpy.array([-1 if n % 2 else +1 for n in numpy.arange(dftsize)])

        else:

            # amplitude "demodulation" in time domain
            self.twiddles_synthesis *= 2 / (1 - numpy.cos(numpy.pi * latency))

        # warmup numba
        dfts = numpy.empty((0, dftsize), complex)
        samples = numpy.empty((0), float)
        buffer = self.buffer
        twiddles = self.twiddles_analysis
        even = self.even
        self.__analyze__(dfts, samples, buffer, twiddles, even)

        # warmup numba
        dfts = numpy.empty((0, dftsize), complex)
        samples = numpy.zeros((0), dtype=float)
        twiddles = self.twiddles_synthesis
        self.__synthesize__(dfts, samples, twiddles)

    def reset(self):
        """
        Reset this SDFT plan to its initial state.
        """

        self.buffer.fill(0)

    def sdft(self, samples):
        """
        Estimate the DFT matrix for the given sample array.

        Parameters
        ----------
        samples : ndarray, list, float
            Array of samples.

        Returns
        -------
        dfts : ndarray
            DFT matrix of shape (samples,frequencies).
        """

        samples = numpy.atleast_1d(samples)

        assert samples.ndim == 1, f'Expected 1D array (samples,), got {samples.shape}!'

        buffer = self.buffer
        twiddles = self.twiddles_analysis
        even = self.even

        dfts = numpy.empty((samples.size, self.size), complex)
        self.__analyze__(dfts, samples, buffer, twiddles, even)
        dfts = self.__convolve__(dfts)
        dfts /= 2

        return dfts

    def isdft(self, dfts):
        """
        Synthesize the sample array from the given DFT matrix.

        Parameters
        ----------
        dfts : ndarray
            DFT matrix of shape (samples,frequencies).

        Returns
        -------
        samples : ndarray
            Array of samples.
        """

        dfts = numpy.atleast_2d(dfts)

        assert dfts.ndim == 2, f'Expected 2D array (samples,frequencies), got {dfts.shape}!'

        twiddles = self.twiddles_synthesis

        samples = numpy.zeros(len(dfts), float)
        self.__synthesize__(dfts, samples, twiddles)
        samples *= 2

        return samples

    def __convolve__(self, x):
        """
        Window the specified DFT matrix.
        """

        x = numpy.atleast_2d(x)

        assert x.ndim == 2, f'Expected 2D array (samples,frequencies), got {x.shape}!'

        M, N = x.shape

        window = str(self.window).lower()

        if window in 'hann':

            y = numpy.hstack((
                numpy.conj(x[:, +1][:, None]),
                x,
                numpy.conj(x[:, -2][:, None])))

            middle = y[:, +1:-1]
            left1  = y[:,   :-2]
            right1 = y[:, +2:  ]

            return (0.5 * middle - 0.25 * (left1 + right1)) / N

        if window in 'hamming':

            y = numpy.hstack((
                numpy.conj(x[:, +1][:, None]),
                x,
                numpy.conj(x[:, -2][:, None])))

            middle = y[:, +1:-1]
            left1  = y[:,   :-2]
            right1 = y[:, +2:  ]

            return (0.54 * middle - 0.23 * (left1 + right1)) / N

        if window in 'blackman':

            y = numpy.hstack((
                numpy.conj(x[:, +2][:, None]),
                numpy.conj(x[:, +1][:, None]),
                x,
                numpy.conj(x[:, -2][:, None]),
                numpy.conj(x[:, -3][:, None])))

            middle = y[:, +2:-2]
            left1  = y[:, +1:-3]
            right1 = y[:, +3:-1]
            left2  = y[:,   :-4]
            right2 = y[:, +4:  ]

            return (0.42 * middle - 0.25 * (left1 + right1) + 0.04 * (left2 + right2)) / N

        return x / N

    @staticmethod
    @numba.jit(nopython=True, fastmath=True)
    def __analyze__(dfts, samples, buffer, twiddles, even):

        if not samples.size:
            return

        first, last = (0.5, 0.5) if even else (0.5, 1.0)

        damping = 1 / twiddles.size

        for i in range(samples.size):

            feedback  = numpy.real(buffer[1:-1]).sum()
            feedback += numpy.real(buffer[0]) * first
            feedback += numpy.real(buffer[-1]) * last
            feedback *= damping

            dfts[i] = (samples[i] - feedback + buffer) * twiddles

            buffer[:] = dfts[i]

    @staticmethod
    @numba.jit(nopython=True, fastmath=True)
    def __synthesize__(dfts, samples, twiddles):

        if not dfts.size:
            return

        for i in range(dfts.shape[0]):

            for j in range(dfts.shape[1]):

                samples[i] += numpy.real(dfts[i, j] * twiddles[j])
