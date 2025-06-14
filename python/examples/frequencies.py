# /// script
# dependencies = ["matplotlib", "numba", "numpy"]
# ///


import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

import matplotlib.pyplot as plot
import numpy as np

from sdft import SDFT


for m in [100, 101]:

    sr = 44100
    w  = None

    sdft  = SDFT(sr, m, w)
    freqs = sdft.frequencies

    d = 1
    t = np.arange(int(d * sr)) / sr

    for n in range(m - m % 2): # skip nyquist

        f = freqs[n]

        print(f'test frequency {f} Hz at bin {n} of {m}')

        x = np.exp(2j * np.pi * f * t).real
        dft = sdft.sdft(x)
        y = sdft.isdft(dft)

        p = np.angle(dft[:, n])
        p = np.unwrap(p)
        p = np.diff(p)
        p = sr * p / (2 * np.pi)

        assert np.allclose(p[p.size//2:], f)
