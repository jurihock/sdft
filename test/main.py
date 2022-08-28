import numpy as np
import os
import sys

from plot import spectrogram, figure, show
from wav import readwav


def main():

    if len(sys.argv) < 5:
        exit(1)

    dftsize = int(sys.argv[1])
    hopsize = int(sys.argv[2])

    ifile = sys.argv[3]
    ofile = sys.argv[4]

    ofiles = {
        'c':   f'{ofile.format("c")}',
        'cpp': f'{ofile.format("cpp")}',
        'py':  f'{ofile.format("py")}'
    }

    x, sr = readwav(ifile)

    y = {
        key: np.fromfile(val, complex).reshape((-1, dftsize))
        for key, val in ofiles.items()
    }

    shapes = list(set([dfts.shape for dfts in y.values()]))
    assert len(shapes) == 1, f'{shapes}'

    assert np.allclose(y['c'], y['cpp']), np.abs(y['c'] - y['cpp']).flatten().max()
    assert np.allclose(y['c'], y['py'], atol=1e-7), np.abs(y['c'] - y['py']).flatten().max()

    figure('c').spectrogram(y['c'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    figure('py').spectrogram(y['py'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    show()


if __name__ == '__main__':

    main()
