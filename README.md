# Sliding Discrete Fourier Transform (SDFT)

According to [[1]](#1) and [[2]](#2).

## Usage

### C

```c
#define SDFT_TD_FLOAT  // time domain data type (float by default)
#define SDFT_FD_DOUBLE // frequency domain data type (double by default)

#include <sdft/sdft.h> // see also src/c folder

size_t n = ...; // number of samples
size_t m = ...; // number of dft bins

float* input = ...; // analysis samples
float* output = ...; // synthesis samples

double complex* dft = ...; // dft matrix of shape (n, m)

sdft_t* sdft = sdft_alloc(m); // create sdft plan

sdft_sdft_n(sdft, n, input, dft); // extract dft matrix from input samples
sdft_isdft_n(sdft, n, dft, output); // synthesize output samples from dft matrix

sdft_free(sdft); // destroy sdft plan
```

### C++

```c++
#include <sdft/sdft.h> // see also src/cpp folder

size_t n = ...; // number of samples
size_t m = ...; // number of dft bins

float* input = ...; // analysis samples
float* output = ...; // synthesis samples

std::complex<double>* dft = ...; // dft matrix of shape (n, m)

SDFT<float, double> sdft(m); // create sdft plan with time/frequency domain data type

sdft.sdft(n, input, dft); // extract dft matrix from input samples
sdft.isdft(n, dft, output); // synthesize output samples from dft matrix
```

## References

1. <span id="1">Krzysztof Duda (2010). Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/5563098</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>
