import numpy

from numpy.lib.stride_tricks import sliding_window_view


def stft(samples, framesize, hopsize):

    frames = sliding_window_view(samples, framesize, writeable=False)[::hopsize]

    N, M = frames.shape

    dfts = numpy.zeros((N, M//2+1), complex)

    w = 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(framesize) / framesize)

    for i, frame in enumerate(frames):

        dfts[i] = numpy.fft.rfft(w * frame, norm='forward')

    return dfts


def istft(dfts, framesize, hopsize):

    N, M = dfts.shape

    samples = numpy.zeros((N * hopsize + framesize), float)

    frames = sliding_window_view(samples, framesize, writeable=True)[::hopsize]

    w = 0.5 - 0.5 * numpy.cos(2 * numpy.pi * numpy.arange(framesize) / framesize)

    w *= hopsize / numpy.sum(w**2)  # force unity gain

    for i, dft in enumerate(dfts):

        frames[i] += w * numpy.fft.irfft(dft, norm='forward')

    return samples
