import numpy

from numpy.lib.stride_tricks import sliding_window_view


class STFT:

    def __init__(self, framesize, hopsize, window='hann', shift=False):

        self.framesize = framesize
        self.hopsize = hopsize
        self.window = window
        self.shift = shift

    def stft(self, samples):

        samples = numpy.atleast_1d(samples)

        assert samples.ndim == 1, f'Expected 1D array (samples,), got {samples.shape}!'

        W = self.weights()

        frames = sliding_window_view(samples, self.framesize, writeable=False)[::self.hopsize]

        dfts = self.fft(frames * W)

        return dfts

    def istft(self, dfts):

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

        if self.shift:

            data = numpy.fft.fftshift(data, axes=-1)

        return numpy.fft.rfft(data, axis=-1, norm='forward')

    def ifft(self, data):

        data = numpy.fft.irfft(data, axis=-1, norm='forward')

        if self.shift:

            return numpy.fft.ifftshift(data, axes=-1)

        return data

    def weights(self):

        N, W = self.framesize, str(self.window).lower()

        if W in 'hann':

            return 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(N) / N)

        if W in 'hamming':

            return 0.54 - 0.46 * numpy.cos(2 * numpy.pi * numpy.arange(N) / N)

        if W in 'blackman':

            return 0.42 - 0.5  * numpy.cos(2 * numpy.pi * numpy.arange(N) / N) \
                        + 0.08 * numpy.cos(4 * numpy.pi * numpy.arange(N) / N)

        return numpy.ones(N)
