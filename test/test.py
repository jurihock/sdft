import numpy as np
import os
import sys


python = os.path.join(os.path.dirname(__file__), '..', 'src', 'python')

sys.path.insert(0, python)


from sdft import SDFT
from wav import read
from plot import spectrogram


def main():

    dftsize = 512
    hopsize = 100
    file = 'test.wav'

    sdft = SDFT(dftsize)

    x, sr = read(file)

    oldsize = x.size
    newsize = (x.size // hopsize) * hopsize

    x = x[:newsize]
    hops = np.split(x, newsize // hopsize)
    y = np.ndarray((len(hops), dftsize), complex)

    for hop, samples in enumerate(hops):

        print(f'{hop+1}/{len(hops)}')

        dfts = sdft.sdft(samples)
        y[hop] = dfts[0]

    spectrogram(y, sr, hopsize, ylim=(500, 15e3), yscale='log').show()


if __name__ == '__main__':

    main()
