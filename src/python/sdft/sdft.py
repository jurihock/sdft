"""
Copyright (c) 2022 Juergen Hock

SPDX-License-Identifier: MIT

Modulated Sliding DFT implementation according to [1] combined with [2].

[1] Krzysztof Duda
    Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform
    IEEE Signal Processing Magazine (2010)
    https://ieeexplore.ieee.org/document/5563098

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

        self.size = dftsize
        self.window = window
        self.latency = latency

        self.offset = 0
        self.delayline = numpy.zeros(dftsize * 2, float)
        self.accumulator = numpy.zeros(dftsize, complex)

        self.twiddles_analysis = numpy.exp(-2j * numpy.pi * numpy.arange(dftsize) / (dftsize * 2))
        self.twiddles_synthesis = numpy.exp(-1j * numpy.pi * numpy.arange(dftsize) * latency)

        if self.latency == 1:

            # circular shift in time domain or multiplication of each dft bin by (-1)**n
            # 1) numpy.array([-1 if n % 2 else +1 for n in numpy.arange(dftsize)])
            # 2) numpy.exp(-1j * numpy.pi * numpy.arange(dftsize) * latency)

            # skip amplitude "demodulation" in time domain
            # since 2 / (1 - numpy.cos(numpy.pi * latency)) is 1

            pass

        else:

            # apply amplitude "demodulation" in time domain
            self.twiddles_synthesis *= 2 / (1 - numpy.cos(numpy.pi * latency))

        dfts = numpy.zeros((0, dftsize), complex)
        samples = numpy.zeros((0), float)
        samples = self.__synthesize__(dfts, samples, self.twiddles_synthesis)
        assert samples.shape == (0, )

    def reset(self):
        """
        Reset this SDFT plan to its initial state.
        """

        self.offset = 0
        self.delayline.fill(0)
        self.accumulator.fill(0)

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

        M = samples.size
        N = self.size

        m = self.offset
        n = numpy.arange(N)

        self.offset += M

        twiddles = self.twiddles_analysis[None, :]
        twiddles = numpy.repeat(twiddles, M + 1, axis=0)
        twiddles[0] **= m
        numpy.cumprod(twiddles, axis=0, out=twiddles)

        delayline = numpy.concatenate((self.delayline, samples))
        numpy.copyto(self.delayline, delayline[-(N * 2):])
        data = samples - delayline[:M]
        data = data[:, None] * twiddles[:-1]

        data[0] += self.accumulator
        numpy.cumsum(data, axis=0, out=data)
        numpy.copyto(self.accumulator, data[-1])
        data *= numpy.conj(twiddles[1:])

        dfts = self.convolve(data)

        return dfts / 2

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

        return self.__synthesize__(dfts, samples, twiddles)

    def convolve(self, x):
        """
        Window the specified DFT matrix.
        """

        window = str(self.window).lower()

        x = numpy.atleast_2d(x)

        assert x.ndim == 2, f'Expected 2D array (samples,frequencies), got {x.shape}!'

        M, N = x.shape

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
    def __synthesize__(dfts, samples, twiddles):

        if not dfts.size:
            return samples

        for i in range(dfts.shape[0]):

            for j in range(dfts.shape[1]):

                samples[i] += numpy.real(dfts[i, j] * twiddles[j])

            samples[i] *= 2

        return samples
