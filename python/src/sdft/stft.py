"""
Copyright (c) 2022 Juergen Hock

SPDX-License-Identifier: MIT

Short-Time Fourier Transform.

Source: https://github.com/jurihock/sdft
"""


import numpy

from numpy.lib.stride_tricks import sliding_window_view


class STFT:
    """
    Short-Time Fourier Transform (STFT).
    """

    def __init__(self, framesize, hopsize, window='hann', shift=False):
        """
        Create a new STFT plan.

        Parameters
        ----------
        framesize : int
            Buffer size in samples.
        hopsize : int
            Hop size in samples.
        window : str, optional
            Analysis window type (boxcar, hann, hamming or blackman).
        shift : bool, optional
            Enable circular shift.
        """

        self.framesize = framesize
        self.hopsize = hopsize
        self.window = window
        self.shift = shift

    def stft(self, samples):
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

        W = self.weights()

        frames = sliding_window_view(samples, self.framesize, writeable=False)[::self.hopsize]

        dfts = self.fft(frames * W)

        return dfts

    def istft(self, dfts):
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

        N, W = dfts.shape[0] * self.hopsize + self.framesize, self.weights()

        W *= self.hopsize / numpy.sum(W**2)  # unity gain

        samples = numpy.zeros((N), float)

        frames0 = sliding_window_view(samples, self.framesize, writeable=True)[::self.hopsize]
        frames1 = self.ifft(dfts) * W

        for i in range(min(len(frames0), len(frames1))):

            frames0[i] += frames1[i]

        return samples

    def fft(self, data):
        """
        Perform forward FFT.
        """

        if self.shift:

            data = numpy.fft.fftshift(data, axes=-1)

        return numpy.fft.rfft(data, axis=-1, norm='forward')

    def ifft(self, data):
        """
        Perform backward FFT.
        """

        data = numpy.fft.irfft(data, axis=-1, norm='forward')

        if self.shift:

            return numpy.fft.ifftshift(data, axes=-1)

        return data

    def weights(self):
        """
        Compute time-domain window coefficients.
        """

        N, W = self.framesize, str(self.window).lower()

        if W in 'hann':

            return 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(N) / N)

        if W in 'hamming':

            return 0.54 - 0.46 * numpy.cos(2 * numpy.pi * numpy.arange(N) / N)

        if W in 'blackman':

            return 0.42 - 0.5  * numpy.cos(2 * numpy.pi * numpy.arange(N) / N) \
                        + 0.08 * numpy.cos(4 * numpy.pi * numpy.arange(N) / N)

        return numpy.ones(N)
