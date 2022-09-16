import numpy

from numpy.lib.stride_tricks import sliding_window_view


def stft(samples, framesize, hopsize, window=None):

    frames = sliding_window_view(samples, framesize, writeable=False)[::hopsize]

    N, M = frames.shape

    dfts = numpy.zeros((N, M//2+1), complex)

    w = weights(framesize, window)

    for i, frame in enumerate(frames):

        dfts[i] = numpy.fft.rfft(w * frame, norm='forward')

    return dfts


def istft(dfts, framesize, hopsize, window=None):

    N, M = dfts.shape

    samples = numpy.zeros((N * hopsize + framesize), float)

    frames = sliding_window_view(samples, framesize, writeable=True)[::hopsize]

    w = weights(framesize, window)
    w *= hopsize / numpy.sum(w**2)

    for i, dft in enumerate(dfts):

        frames[i] += w * numpy.fft.irfft(dft, norm='forward')

    return samples


def weights(size, window=None):

    window = str(window).lower()

    if window in 'hann':

        return 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

    if window in 'hamming':

        return 0.54 - 0.46 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

    if window in 'blackman':

        return 0.42 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size) \
                    + 0.08 * numpy.cos(4 * numpy.pi * numpy.arange(size) / size)

    return numpy.ones(size)
