import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

import matplotlib.pyplot as plot
import numpy as np

from sdft import STFT

asymmetric = False
overlap = 1

winsize = 128
hopsize = 128 // overlap
dftsize = 128+1 if asymmetric else None

stft = STFT(winsize, hopsize, dftsize)

print('A', stft.analysis_window_size, 'S', stft.synthesis_window_size)

x = np.ones(winsize * 10)
y = stft.istft(stft.stft(x))

h = stft.symmetric_window(winsize)
h = np.tile(h, len(y) // winsize)
h = (h**2) * np.max(y)

plot.title('asymmetric' if asymmetric else 'symmetric')
plot.plot(y, color='b', alpha=0.8, label='wola')
plot.plot(h, color='r', alpha=0.8, label='hann')
plot.legend()
plot.show()
