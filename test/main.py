import numpy as np
import os
import sys

from plot import spectrogram, figure, show
from wav import readwav


def main():

    dftsize = 512
    hopsize = 1000

    ifile = 'test.wav'
    ofile = 'test.*.dfts'

    ofiles = {
        'c': f'../build/{ofile.replace("*", "c")}',
        'cpp': f'../build/{ofile.replace("*", "cpp")}',
        'py': f'{ofile.replace("*", "py")}'
    }

    x, sr = readwav(ifile)

    y = {
        key: np.fromfile(val, complex).reshape((-1, dftsize))
        for key, val in ofiles.items()
    }

    shapes = list(set([dfts.shape for dfts in y.values()]))
    assert len(shapes) == 1, f'{shapes}'

    figure('c').spectrogram(y['c'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    figure('cpp').spectrogram(y['cpp'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    figure('py').spectrogram(y['py'], sr, hopsize, ylim=(500, 15e3), yscale='log')
    show()


if __name__ == '__main__':

    main()
