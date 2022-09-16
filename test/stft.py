import numpy

from numpy.lib.stride_tricks import sliding_window_view


def stft(samples, framesize, hopsize, window='hann', shift=False):

    frames = sliding_window_view(samples, framesize, writeable=False)[::hopsize]

    N, M = frames.shape

    dfts = numpy.zeros((N, M//2+1), complex)

    w = weights(framesize, window)

    for i, frame in enumerate(frames):

        dfts[i] = fft(w * frame, shift)

    return dfts


def istft(dfts, framesize, hopsize, window='hann', shift=False):

    N, M = dfts.shape

    samples = numpy.zeros((N * hopsize + framesize), float)

    frames = sliding_window_view(samples, framesize, writeable=True)[::hopsize]

    w = weights(framesize, window)
    w *= hopsize / numpy.sum(w**2)

    for i, dft in enumerate(dfts):

        frames[i] += w * ifft(dft, shift)

    return samples


def fft(data, shift=False):

    if shift:

        return numpy.fft.rfft(numpy.fft.fftshift(data), norm='forward')

    else:

        return numpy.fft.rfft(data, norm='forward')


def ifft(data, shift=False):

    if shift:

        return numpy.fft.fftshift(numpy.fft.irfft(data, norm='forward'))

    else:

        return numpy.fft.irfft(data, norm='forward')


def weights(size, window='hann'):

    window = str(window).lower()

    if window in 'hann':

        return 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

    if window in 'hamming':

        return 0.54 - 0.46 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size)

    if window in 'blackman':

        return 0.42 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(size) / size) \
                    + 0.08 * numpy.cos(4 * numpy.pi * numpy.arange(size) / size)

    return numpy.ones(size)
