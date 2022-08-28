import numpy as np
import os
import sys


sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..', 'src', 'python'))


from sdft import SDFT
from wav import readwav
from plot import spectrogram


def main():

    dftsize = 512
    hopsize = 1000

    ifile = 'test.wav'
    ofile = 'test.py.dfts'

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
