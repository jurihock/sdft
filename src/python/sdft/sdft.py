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


import numpy


class SDFT:
    """
    Sliding Discrete Fourier Transform (SDFT).
    """

    def __init__(self, dftsize, latency=1):
        """
        Create a new SDFT plan.

        Parameters
        ----------
        dftsize : int
            Desired number of DFT bins.
        latency : float
            Synthesis latency factor between 0 and 1.
            The default value 1 corresponds to the highest latency and best possible SNR.
            A smaller value decreases both latency and SNR, but also increases the workload.
        """

        self.size = dftsize
        self.latency = latency

        self.offset = 0
        self.delayline = numpy.zeros(dftsize * 2, float)
        self.accumulator = numpy.zeros(dftsize, complex)

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

        m = numpy.arange(self.offset, self.offset + M + 1)[:, None]
        n = numpy.arange(N)

        self.offset += M

        twiddles = numpy.exp(-2j * numpy.pi * m * n / (N * 2))

        delayline = numpy.concatenate((self.delayline, samples))
        numpy.copyto(self.delayline, delayline[-(N * 2):])
        data = samples - delayline[:M]
        data = data[:, None] * twiddles[:-1]

        data[0] += self.accumulator
        numpy.add.accumulate(data, axis=0, out=data)
        numpy.copyto(self.accumulator, data[-1])
        data *= numpy.conj(twiddles[1:])

        dfts = self.window(data)

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

        M, N = dfts.shape

        twiddles = numpy.array([-1 if n % 2 else +1 for n in numpy.arange(N)]) \
                   if self.latency == 1 else \
                   numpy.exp(-1j * numpy.pi * self.latency * numpy.arange(N))

        weight = 2 / (1 - numpy.cos(numpy.pi * self.latency))

        samples = numpy.sum(numpy.real(dfts * twiddles * weight), axis=-1)

        samples *= 2

        return samples

    def window(self, x):
        """
        Apply Hann window to the specified DFT matrix.
        """

        x = numpy.atleast_2d(x)

        assert x.ndim == 2, f'Expected 2D array (samples,frequencies), got {x.shape}!'

        M, N = x.shape

        y = numpy.hstack((
            numpy.conj(x[:,1][:,None]),
            x,
            numpy.conj(x[:,-2][:,None])))

        left = y[:, :-2]
        right = y[:, +2:]
        middle = y[:, +1:-1]

        return ((middle + middle) - (left + right)) / (N * 4)
