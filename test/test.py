import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..', 'src', 'python'))

from sdft import SDFT
from wav import readwav, writewav


def main():

    if len(sys.argv) < 6:
        exit(1)

    dftsize = int(sys.argv[1])
    hopsize = int(sys.argv[2])
    srcfile = sys.argv[3]
    wavfile = sys.argv[4]
    dftfile = sys.argv[5]

    sdft = SDFT(dftsize)

    x, sr = readwav(srcfile)
    size = x.size

    print(f'PY\t{srcfile} {size} {sr}');
    size = (x.size // hopsize) * hopsize

    x = x[:size]
    y = np.ndarray((0), float)
    hops = np.split(x, size // hopsize)
    dfts = np.ndarray((len(hops), dftsize), complex)

    progress = 0

    for hop, samples in enumerate(hops):

        percent = (hop + 1) / len(hops)

        if int(percent * 10) != progress:
            progress = int(percent * 10)
            print(f'{progress * 10}%')

        tmp = sdft.sdft(samples)
        dfts[hop] = tmp[0]
        y = np.concatenate((y, sdft.isdft(tmp)))

    writewav(wavfile, y, sr)
    dfts.tofile(dftfile)


if __name__ == '__main__':

    main()
