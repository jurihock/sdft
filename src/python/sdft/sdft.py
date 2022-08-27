import numpy


def window(x):

    M, N = x.shape

    left = x[:, :-2]
    right = x[:, +2:]
    middle = x[:, +1:-1]

    y = ((middle + middle) - (left + right)) / (N * 4)

    y = numpy.pad(y, ((0, 0), (1, 1)))

    return y


def sdft(samples, dftsize):

    M = samples.size
    N = dftsize

    m = numpy.arange(M + 1)[:, None]
    n = numpy.arange(N)

    twiddles = numpy.exp(-2j * numpy.pi * m * n / (N * 2))

    delayline = numpy.resize(numpy.pad(samples, (N * 2, 0)), samples.shape)

    deltas = samples - delayline

    data = deltas[:, None] * twiddles[:-1]

    numpy.add.accumulate(data, axis=0, out=data)

    data *= numpy.conj(twiddles[1:])

    dfts = window(data)

    return dfts


def isdft(dfts, latency=1):

    M, N = dfts.shape

    twiddles = numpy.array([-1 if n % 2 else +1 for n in numpy.arange(N)]) \
               if latency == 1 else \
               numpy.exp(-1j * numpy.pi * latency * numpy.arange(N))

    weight = 2 / (1 - numpy.cos(numpy.pi * latency))

    samples = numpy.sum(numpy.real(dfts * twiddles * weight), axis=-1)

    return samples
