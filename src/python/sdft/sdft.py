import numpy


class SDFT:

    def __init__(self, dftsize, latency=1):

        self.dftsize = dftsize
        self.latency = latency

        self.offset = 0
        self.delayline = numpy.zeros(dftsize * 2, float)
        self.accumulator = numpy.zeros(dftsize, complex)

    def sdft(self, samples):

        samples = numpy.atleast_1d(samples)

        assert samples.ndim == 1, f'Expected 1D array (samples,), got {samples.shape}!'

        M = samples.size
        N = self.dftsize

        m = numpy.arange(self.offset, self.offset + M + 1)[:, None]
        n = numpy.arange(N)

        self.offset += M

        twiddles = numpy.exp(-2j * numpy.pi * m * n / (N * 2))

        delayline = numpy.concatenate((self.delayline, samples))
        deltas = samples - delayline[:M]
        numpy.copyto(self.delayline, delayline[-(N * 2):])

        data = deltas[:, None] * twiddles[:-1]

        data[0] += self.accumulator
        numpy.add.accumulate(data, axis=0, out=data)
        numpy.copyto(self.accumulator, data[-1])

        data *= numpy.conj(twiddles[1:])
        dfts = self.window(data)

        return dfts

    def isdft(self, dfts):

        dfts = numpy.atleast_2d(dfts)

        assert dfts.ndim == 2, f'Expected 2D array (samples,frequencies), got {dfts.shape}!'

        M, N = dfts.shape

        twiddles = numpy.array([-1 if n % 2 else +1 for n in numpy.arange(N)]) \
                   if self.latency == 1 else \
                   numpy.exp(-1j * numpy.pi * self.latency * numpy.arange(N))

        weight = 2 / (1 - numpy.cos(numpy.pi * self.latency))

        samples = numpy.sum(numpy.real(dfts * twiddles * weight), axis=-1)

        return samples

    def window(self, x):

        x = numpy.atleast_2d(x)

        assert x.ndim == 2, f'Expected 2D array (samples,frequencies), got {x.shape}!'

        M, N = x.shape

        left = x[:, :-2]
        right = x[:, +2:]
        middle = x[:, +1:-1]

        y = ((middle + middle) - (left + right)) / (N * 4)

        y = numpy.pad(y, ((0, 0), (1, 1)))

        return y
