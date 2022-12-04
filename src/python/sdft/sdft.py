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

    def __init__(self, dftsize, window='hann', latency=None):
        """
        Create a new SDFT plan.

        Parameters
        ----------
        dftsize : int
            Desired number of DFT bins.
        window : str, optional
            Analysis window type (boxcar, hann, hamming or blackman).
        latency : float, optional
            Number of samples to perform time-shifting and thus affect the output latency.
            The specified value can be a positive or a negative number.
            The default value None disables latency compensation.
            The special value 0 enables zero-latency output.
        """

        self.size = dftsize
        self.window = window
        self.latency = latency

        self.offset = 0
        self.delayline = numpy.zeros(dftsize * 2, float)
        self.accumulator = numpy.zeros(dftsize, complex)

        self.omega = -2j * numpy.pi / (dftsize * 2)
        self.twiddles = numpy.exp(self.omega * numpy.arange(dftsize))
        self.timeshift = self.latency if self.latency else -(dftsize - 1) if self.latency == 0 else None
        self.timeshift = numpy.exp(self.omega * numpy.arange(dftsize) * self.timeshift) if self.timeshift else None

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

        twiddles = self.twiddles[None, :]
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

        if self.timeshift is not None:
            data *= self.timeshift[None, :]

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

        M, N = dfts.shape

        twiddles = numpy.array([-1 if n % 2 else +1 for n in numpy.arange(N)])
        samples = numpy.sum(numpy.real(dfts * twiddles), axis=-1)

        return samples * 2

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
                numpy.conj(x[:,+1][:,None]),
                x,
                numpy.conj(x[:,-2][:,None])))

            middle = y[:, +1:-1]
            left1  = y[:,   :-2]
            right1 = y[:, +2:  ]

            return (0.5 * middle - 0.25 * (left1 + right1)) / N

        if window in 'hamming':

            y = numpy.hstack((
                numpy.conj(x[:,+1][:,None]),
                x,
                numpy.conj(x[:,-2][:,None])))

            middle = y[:, +1:-1]
            left1  = y[:,   :-2]
            right1 = y[:, +2:  ]

            return (0.54 * middle - 0.23 * (left1 + right1)) / N

        if window in 'blackman':

            y = numpy.hstack((
                numpy.conj(x[:,+2][:,None]),
                numpy.conj(x[:,+1][:,None]),
                x,
                numpy.conj(x[:,-2][:,None]),
                numpy.conj(x[:,-3][:,None])))

            middle = y[:, +2:-2]
            left1  = y[:, +1:-3]
            right1 = y[:, +3:-1]
            left2  = y[:,   :-4]
            right2 = y[:, +4:  ]

            return (0.42 * middle - 0.25 * (left1 + right1) + 0.04 * (left2 + right2)) / N

        return x / N
