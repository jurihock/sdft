import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..', 'src', 'python'))

from sdft import SDFT
from wav import readwav


def main():

    if len(sys.argv) < 5:
        exit(1)

    dftsize = int(sys.argv[1])
    hopsize = int(sys.argv[2])

    ifile = sys.argv[3]
    ofile = sys.argv[4]

    sdft = SDFT(dftsize)

    x, sr = readwav(ifile)
    size = x.size

    print(f'PY\t{ifile} {size} {sr}');
    size = (x.size // hopsize) * hopsize

    x = x[:size]
    hops = np.split(x, size // hopsize)

    y = np.ndarray((len(hops), dftsize), complex)

    for hop, samples in enumerate(hops):

        print(f'{hop+1}/{len(hops)}')
        y[hop] = sdft.sdft(samples)[0]

    y.tofile(ofile)


if __name__ == '__main__':

    main()
