import numpy as np
import os
import sys

from plot import spectrogram, figure, show
from wav import readwav


def main():

    if len(sys.argv) < 6:
        exit(1)

    dftsize = int(sys.argv[1])
    hopsize = int(sys.argv[2])

    srcfile = sys.argv[3]
    wavfile = sys.argv[4]
    dftfile = sys.argv[5]

    wavfiles = {
        'c':   f'{wavfile.format("c")}',
        'cpp': f'{wavfile.format("cpp")}',
        'py':  f'{wavfile.format("py")}'
    }

    dftfiles = {
        'c':   f'{dftfile.format("c")}',
        'cpp': f'{dftfile.format("cpp")}',
        'py':  f'{dftfile.format("py")}'
    }

    x, sr = readwav(srcfile)

    y = {
        key: readwav(val)[0]
        for key, val in wavfiles.items()
    }

    dfts = {
        key: np.fromfile(val, complex).reshape((-1, dftsize))
        for key, val in dftfiles.items()
    }

    # check wavs

    shapes = list(set([_.shape for _ in y.values()]))
    assert len(shapes) == 1, f'{shapes}'

    assert np.allclose(y['c'], y['cpp']), np.abs(y['c'] - y['cpp']).flatten().max()
    # TODO assert np.allclose(y['c'], y['py']), np.abs(y['c'] - y['py']).flatten().max()

    # check dfts

    shapes = list(set([_.shape for _ in dfts.values()]))
    assert len(shapes) == 1, f'{shapes}'

    assert np.allclose(dfts['c'], dfts['cpp']), np.abs(dfts['c'] - dfts['cpp']).flatten().max()
    assert np.allclose(dfts['c'], dfts['py'], atol=1e-7), np.abs(dfts['c'] - dfts['py']).flatten().max()

    # plot spectrograms

    figure('c').spectrogram(dfts['c'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    figure('py').spectrogram(dfts['py'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    show()


if __name__ == '__main__':

    main()
