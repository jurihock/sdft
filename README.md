# Sliding Discrete Fourier Transform (SDFT)

![language](https://img.shields.io/badge/languages-C%2FC%2B%2B%20Python-blue)
![license](https://img.shields.io/github/license/jurihock/sdft?color=green)
![pypi](https://img.shields.io/pypi/v/sdft?color=gold)

Forward and inverse Sliding DFT according to [[1]](#1) and [[2]](#2) with following features:

- Arbitrary number of DFT bins
- Hann window at analysis
- Single or multiple sample processing at once
- Endless forward/inverse transform
- Optional synthesis latency control parameter

The [Sliding Discrete Fourier Transform (SDFT)](https://en.wikipedia.org/wiki/Sliding_DFT) is a recursive approach to compute the Fourier transform sample by sample. In this particular case it's more efficient than the FFT based [Short Time Fourier Transform (STFT)](https://en.wikipedia.org/wiki/Short-time_Fourier_transform) approach with one sample hops. On the other side, the SDFT is still known to suffer from accumulated errors and potential instabilities.

This implementation features the *modulated* SDFT algorithm, which is guaranteed to be stable while being accurate. It takes real valued samples and estimates the corresponding half size complex valued DFT vector for each of them. The length of the estimated DFT vector is not limited to the power of two. The eventually altered DFT vector can also be used to synthesize an output sample.

Compared to STFT, the algorithmic synthesis latency of SDFT is lower and can additionally be reduced at the expense of signal to noise ratio. This is particularly advantageous for real-time applications, e.g. digital audio signal processing.

## Usage

### C

```c
#define SDFT_TD_FLOAT  // time domain data type (float by default)
#define SDFT_FD_DOUBLE // frequency domain data type (double by default)

#include <sdft/sdft.h> // see also src/c folder

size_t n = ...; // number of samples
size_t m = ...; // number of dft bins

float* x = ...; // analysis samples of shape (n)
float* y = ...; // synthesis samples of shape (n)

double complex* dft = ...; // dft matrix of shape (n, m)

sdft_t* sdft = sdft_alloc(m); // create sdft plan

sdft_sdft_n(sdft, n, x, dft); // extract dft matrix from input samples
sdft_isdft_n(sdft, n, dft, y); // synthesize output samples from dft matrix

sdft_free(sdft); // destroy sdft plan
```

Due to incomplete C complex math support in MSVC, optionally use following universal typedefs:

* `sdft_float_t` instead of `float`
* `sdft_double_complex_t` instead of `double complex`

or even better the corresponding generic typedefs:

* `sdft_td_t`
* `sdft_fdx_t`

In both cases, the underlying data type results from the `SDFT_TD_*` and `SDFT_FD_*` definitions.

### C++

```c++
#include <sdft/sdft.h> // see also src/cpp folder

size_t n = ...; // number of samples
size_t m = ...; // number of dft bins

float* x = ...; // analysis samples of shape (n)
float* y = ...; // synthesis samples of shape (n)

std::complex<double>* dft = ...; // dft matrix of shape (n, m)

SDFT<float> sdft(m); // create sdft plan with time/frequency domain data type

sdft.sdft(n, x, dft); // extract dft matrix from input samples
sdft.isdft(n, dft, y); // synthesize output samples from dft matrix
```

Optionally specify a different frequency domain data type like so `SDFT<float, double>`.

### Python

```python
from sdft import SDFT # see also src/python folder

n = ... # number of samples
m = ... # number of dft bins

x = ... # analysis samples of shape (n)

sdft = SDFT(m) # create sdft plan

dft = sdft.sdft(x) # extract dft matrix from input samples
y = sdft.isdft(dft) # synthesize output samples from dft matrix
```

Feel free to obtain the current version from the PyPI `pip install sdft`.

## Test spectrogram

Below you can see two spectrograms of the same file `test.wav` computed by SDFT and STFT with identical spectral resolution, window function and hop size. Do you see any significant difference between them?

| SDFT | STFT |
| ---- | ---- |
| ![SDFT](https://github.com/jurihock/sdft/blob/main/test/sdft.png) | ![STFT](https://github.com/jurihock/sdft/blob/main/test/stft.png) |

Well, the results are very similar, which is to be considered as the proof of concept...

## Implementation details

- The C/C++ implementation corresponds to figure 4 in [[1]](#1).
- The Python implementation corresponds to figure 3b in [[1]](#1).

## References

1. <span id="1">Krzysztof Duda (2010). Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/5563098</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>

## License

This *sdft* implementation is licensed under the terms of the MIT license.
For details please refer to the accompanying [LICENSE](LICENSE) file distributed with it.
