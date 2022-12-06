import numpy as np
import os
import sys

from plot import spectrogram, figure, show
from stft import STFT
from wav import readwav


def main():

    if len(sys.argv) < 8:
        exit(1)

    dftsize = int(sys.argv[1])
    hopsize = int(sys.argv[2])
    window = sys.argv[3]
    latency = float(sys.argv[4])
    srcfile = sys.argv[5]
    wavfile = sys.argv[6]
    dftfile = sys.argv[7]

    framesize = dftsize * 2

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

    # emulate real-time delay of one frame

    x = np.roll(x, framesize)
    x[:framesize] = 0

    # compute stft reference spectrogram

    stft = STFT(framesize, hopsize, window).stft(x)

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

    figure('stft').spectrogram(stft, sr, hopsize, xlim=(0, 7.9), ylim=(500, 15e3), yscale='log').tight()
    figure('c').spectrogram(dfts['c'], sr, hopsize, xlim=(0, 7.9), ylim=(500, 15e3), yscale='log').tight()
    figure('py').spectrogram(dfts['py'], sr, hopsize, xlim=(0, 7.9), ylim=(500, 15e3), yscale='log').tight()
    show()


if __name__ == '__main__':

    main()
