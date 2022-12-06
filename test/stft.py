import numpy

from numpy.lib.stride_tricks import sliding_window_view


class STFT:

    def __init__(self, framesize, hopsize, window='hann', shift=False):

        self.framesize = framesize
        self.dftsize = numpy.fft.rfftfreq(framesize).size
        self.hopsize = hopsize
        self.window = window
        self.shift = shift

    def stft(self, samples):

        samples = numpy.atleast_1d(samples)

        assert samples.ndim == 1, f'Expected 1D array (samples,), got {samples.shape}!'

        frames = sliding_window_view(samples, self.framesize, writeable=False)[::self.hopsize]

        M, N, W = frames.shape[0], self.dftsize, self.weights()

        dfts = self.fft(frames * W)

        return dfts

    def istft(self, dfts):

        dfts = numpy.atleast_2d(dfts)

        assert dfts.ndim == 2, f'Expected 2D array (samples,frequencies), got {dfts.shape}!'

        M, W = dfts.shape[0] * self.hopsize + self.framesize, self.weights()

        W *= self.hopsize / numpy.sum(W**2)

        samples = numpy.zeros((M), float)

        frames = sliding_window_view(samples, self.framesize, writeable=True)[::self.hopsize]
        frames += self.ifft(dft) * W

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

        size = self.framesize
        window = str(self.window).lower()

        if window in 'hann':

            return 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

        if window in 'hamming':

            return 0.54 - 0.46 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

        if window in 'blackman':

            return 0.42 - 0.5  * numpy.cos(2 * numpy.pi * numpy.arange(size) / size) \
                        + 0.08 * numpy.cos(4 * numpy.pi * numpy.arange(size) / size)

        return numpy.ones(size)
