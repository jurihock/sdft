import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

import matplotlib.pyplot as plot
import numpy as np

from sdft import SDFT


def signal(n):

    import scipy.stats
    return scipy.stats.truncnorm(-1, +1, loc=0, scale=1).rvs(n)


def main():

    # SUBJECT: modify l and m variables to see what happens


    # 1) generate input signal

    sr = 44100            # sample rate in Hz
    n = int(100e-3 * sr)  # number of samples
    m = 1000              # number of dft bins
    l = 1                 # latency reduction factor in range (0..1]
                          # (which is 1 by default)

    c = int((m - 1) * l)  # latency in samples
    n += c                # compensate latency

    t = np.arange(n) / sr  # timestamps in seconds

    x = signal(n)  # sample vector of shape (n)


    # 2) synthesize output signal

    sdft = SDFT(m, latency=l)  # create sdft plan

    dft = sdft.sdft(x)  # dft matrix of shape (n, m)

    y = sdft.isdft(dft)  # sample vector of shape (n)


    # 3) plot x and y signals

    y = y[c:]
    x = x[:y.size]
    t = t[:x.size]

    e = y - x  # error signal or "noise"

    snr = np.mean(x**2) / np.mean(e**2)  # signal to noise (error) ratio
    snr = int(10 * np.log10(snr))  # convert to decibel

    l = int(1e3 * c / sr)  # convert samples to ms
    t *= 1e3  # convert s to ms

    plot.plot(t, x, label='x', alpha=0.9)
    plot.plot(t, y, label='y', alpha=0.9)

    plot.title(f'Latency {l} ms, SNR {snr} dB (Signal / Error)')
    plot.legend()
    plot.xlabel('ms')

    plot.ylim(-1.1, +1.1)

    plot.show()


if __name__ == '__main__':

    main()
