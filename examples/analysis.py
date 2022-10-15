from sdft import SDFT

import matplotlib.pyplot as plot
import numpy as np


def phase(t, f):

    dt = np.diff(t, prepend=0)
    dp = 2 * np.pi * f * dt
    return np.cumsum(dp)


def main():

    # 1) generate input signal

    sr = 44100  # sample rate in Hz
    n = 1 * sr  # number of samples
    m = 1000    # number of dft bins

    t = np.arange(n) / sr  # timestamps in seconds

    f = 1000  # single frequency
    # f = np.linspace(0, 1000, n)  # linear chirp
    # f = np.sin(np.pi * t * t[-1]) * 1000 + 1000  # frequency wave

    x = np.sin(phase(t, f))  # sample vector of shape (n)


    # 2) estimate output dft

    sdft = SDFT(m)  # create sdft plan

    dft = sdft.sdft(x)  # dft matrix of shape (n, m)


    # 3) plot spectrogram

    with np.errstate(all='ignore'):
        db = 20 * np.log10(np.abs(dft))

    roi = (0, n / sr, 0, sr / 2)
    args = dict(extent=roi, origin='lower', aspect='auto', cmap='inferno', interpolation='nearest')

    plot.imshow(db.T, **args)
    cbar = plot.colorbar()

    plot.xlabel('s')
    plot.ylabel('Hz')
    cbar.set_label('dB')

    plot.ylim(0, 5000)
    plot.clim(-120, 0)

    plot.show()


if __name__ == '__main__':

    main()
