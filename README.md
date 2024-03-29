# Sliding Discrete Fourier Transform (SDFT)

![language](https://img.shields.io/badge/languages-C%2FC%2B%2B%20Python-blue)
![license](https://img.shields.io/github/license/jurihock/sdft?color=green)
![pypi](https://img.shields.io/pypi/v/sdft?color=gold)

Forward and inverse Sliding DFT according to [[1]](#1) and [[2]](#2) with following features:

- Arbitrary number of DFT bins
- Built-in analysis window functions Boxcar, Hann (default), Hamming and Blackman
- Customizable time and frequency domain data type in C/C++
- Endless single or multiple sample processing at once
- Optional synthesis latency control parameter
- Real-time low latency<sup>*</sup> analysis and synthesis capability

<sup><sub>*) compared to STFT latency</sub></sup>

The [Sliding Discrete Fourier Transform (SDFT)](https://en.wikipedia.org/wiki/Sliding_DFT) is a recursive approach to compute the Fourier transform sample by sample. In this particular case it's more efficient than the FFT based [Short Time Fourier Transform (STFT)](https://en.wikipedia.org/wiki/Short-time_Fourier_transform) approach with one sample hops. On the other side, the SDFT is still known to suffer from accumulated errors and potential instabilities.

This implementation features the *modulated* SDFT algorithm, which is guaranteed to be stable while being accurate. It takes real valued samples and estimates the corresponding half size complex valued DFT vector for each of them. The DFT vector of a particular input sample contains linearly spaced frequency bins from 0 Hz up to the [Nyquist](https://en.wikipedia.org/wiki/Nyquist_frequency) frequency. Therefore, the spectral resolution depends on the input sample rate and the specified DFT vector length. The length of the estimated DFT vector, as specified in the SDFT constructor, is not limited to the power of two. The eventually altered DFT vector can also be used to synthesize an output sample.

Compared to STFT, the algorithmic synthesis latency of SDFT is lower and can additionally be reduced at the expense of signal to noise ratio. Spectral data processing coupled with reduced latency is especially useful for real-time applications, e.g. digital audio signal processing.

## Basic usage

For detailed usage, please have a look at the provided "analysis" examples.

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

<details>
<summary><strong>MSVC</strong></summary>
<p/>

Due to incomplete [C complex math support](https://docs.microsoft.com/cpp/c-runtime-library/complex-math-support) in MSVC, optionally use following universal typedefs:

* `sdft_float_t` instead of `float`
* `sdft_double_complex_t` instead of `double complex`

or even better the corresponding generic typedefs:

* `sdft_td_t`
* `sdft_fdx_t`

In both cases, the underlying data type results from the `SDFT_TD_*` and `SDFT_FD_*` definitions.

</details>

<details>
<summary><strong>No complex.h? No problem...</strong></summary>
<p/>

Just define `SDFT_NO_COMPLEX_H` to prevent `complex.h` from being included and internally enable compatible complex number representation instead:

```c
typedef struct { sdft_fd_t r, i; } sdft_fdx_t;
```

</details>

### C++

```c++
#include <sdft/sdft.h> // see also src/cpp folder

size_t n = ...; // number of samples
size_t m = ...; // number of dft bins

float* x = ...; // analysis samples of shape (n)
float* y = ...; // synthesis samples of shape (n)

std::complex<double>* dft = ...; // dft matrix of shape (n, m)

SDFT<float, double> sdft(m); // create sdft plan for custom time and frequency domain data types

sdft.sdft(n, x, dft); // extract dft matrix from input samples
sdft.isdft(n, dft, y); // synthesize output samples from dft matrix
```

The time domain data type defaults to `float` and the frequency domain data type to `double`.

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

Feel free to obtain current version from [PyPI](https://pypi.org/project/sdft) by executing `pip install sdft`.

## Test spectrogram

Below you can see two spectrograms of the same audio file `test.wav` computed by SDFT and STFT with identical spectral resolution, window function and hop size. Do you see any significant differences between them?

| SDFT | STFT |
| ---- | ---- |
| ![SDFT](https://github.com/jurihock/sdft/raw/main/test/sdft.png) | ![STFT](https://github.com/jurihock/sdft/raw/main/test/stft.png) |

Well, the results are very similar, which is to be considered as the proof of concept...

## See also

If you're interested in Sliding DFT with *logarithmic* frequency resolution, don't forget to browse my [jurihock/qdft](https://github.com/jurihock/qdft) project!

## References

1. <span id="1">Krzysztof Duda (2010). Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/5563098</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>

## License

[github.com/jurihock/sdft](https://github.com/jurihock/sdft) is licensed under the terms of the MIT license.
For details please refer to the accompanying [LICENSE](https://github.com/jurihock/sdft/raw/main/LICENSE) file distributed with it.
