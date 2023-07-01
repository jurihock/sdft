/**
 * Copyright (c) 2022 Juergen Hock
 *
 * SPDX-License-Identifier: MIT
 *
 * Modulated Sliding DFT implementation according to [1] combined with [2].
 *
 * [1] Krzysztof Duda
 *     Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform
 *     IEEE Signal Processing Magazine (2010)
 *     https://ieeexplore.ieee.org/document/5563098
 *
 * [2] Russell Bradford and Richard Dobson and John ffitch
 *     Sliding is Smoother than Jumping
 *     International Computer Music Conference (2005)
 *     http://hdl.handle.net/2027/spo.bbp2372.2005.086
 *
 * Source: https://github.com/jurihock/sdft
 **/

/**
 * Use one of following defines to specify the desired
 * time domain (TD) and frequency domain (FD) data type.
 *
 * #define SDFT_TD_FLOAT        // default
 * #define SDFT_TD_DOUBLE
 * #define SDFT_TD_LONG_DOUBLE
 *
 * #define SDFT_FD_FLOAT
 * #define SDFT_FD_DOUBLE       // default and recommended
 * #define SDFT_FD_LONG_DOUBLE
 *
 * The specified data type appears as typedef
 * sdft_td_t, sdft_fd_t and sdft_fdx_t for complex numbers.
 *
 * Optionally define SDFT_NO_COMPLEX_H to not include the complex.h header.
 **/

/**
 * List of common functions to start with:
 *
 * sdft_alloc[_custom]
 * sdft_free
 *
 * sdft_reset
 * sdft_size
 * sdft_latency
 *
 * sdft_sdft
 * sdft_sdft_n
 * sdft_sdft_nd
 *
 * sdft_isdft
 * sdft_isdft_n
 * sdft_isdft_nd
 *
 * Functions sdft_etc_* are reserved for internal use.
 **/

#pragma once

#include <math.h>
#include <stdlib.h>
#include <string.h>

#if !defined(SDFT_NO_COMPLEX_H)
#include <complex.h>
#endif

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(_MSC_VER)
#define SDFT_MSVC
#endif

typedef size_t sdft_size_t;

typedef float sdft_float_t;
typedef double sdft_double_t;
typedef long double sdft_long_double_t;

#if defined(SDFT_NO_COMPLEX_H)
  struct sdft_float_complex { float r, i; };
  struct sdft_double_complex { double r, i; };
  struct sdft_long_double_complex { long double r, i; };
  typedef struct sdft_float_complex sdft_float_complex_t;
  typedef struct sdft_double_complex sdft_double_complex_t;
  typedef struct sdft_long_double_complex sdft_long_double_complex_t;
#elif defined(SDFT_MSVC)
  typedef _Fcomplex sdft_float_complex_t;
  typedef _Dcomplex sdft_double_complex_t;
  typedef _Lcomplex sdft_long_double_complex_t;
#else
  typedef float complex sdft_float_complex_t;
  typedef double complex sdft_double_complex_t;
  typedef long double complex sdft_long_double_complex_t;
#endif

#if defined(SDFT_TD_FLOAT)
  typedef sdft_float_t sdft_td_t;
#elif defined(SDFT_TD_DOUBLE)
  typedef sdft_double_t sdft_td_t;
#elif defined(SDFT_TD_LONG_DOUBLE)
  typedef sdft_long_double_t sdft_td_t;
#else
  #define SDFT_TD_FLOAT
  typedef sdft_float_t sdft_td_t;
#endif

#if defined(SDFT_FD_FLOAT)
  typedef sdft_float_t sdft_fd_t;
  typedef sdft_float_complex_t sdft_fdx_t;
#elif defined(SDFT_FD_DOUBLE)
  typedef sdft_double_t sdft_fd_t;
  typedef sdft_double_complex_t sdft_fdx_t;
#elif defined(SDFT_FD_LONG_DOUBLE)
  typedef sdft_long_double_t sdft_fd_t;
  typedef sdft_long_double_complex_t sdft_fdx_t;
#else
  #define SDFT_FD_DOUBLE
  typedef sdft_double_t sdft_fd_t;
  typedef sdft_double_complex_t sdft_fdx_t;
#endif

enum sdft_window
{
  sdft_window_boxcar,
  sdft_window_hann,
  sdft_window_hamming,
  sdft_window_blackman
};

typedef enum sdft_window sdft_window_t;

struct sdft_plan_roi
{
  sdft_size_t first;
  sdft_size_t second;
};

typedef struct sdft_plan_roi sdft_roi_t;

struct sdft_plan_analysis
{
  sdft_window_t window;

  sdft_fd_t weight;
  sdft_roi_t roi;
  sdft_fdx_t* twiddles;

  sdft_size_t cursor;
  sdft_size_t maxcursor;
  sdft_td_t* input;

  sdft_fdx_t* accoutput;
  sdft_fdx_t* auxoutput;
  sdft_fdx_t* fiddles;
};

typedef struct sdft_plan_analysis sdft_analysis_t;

struct sdft_plan_synthesis
{
  sdft_double_t latency;

  sdft_fd_t weight;
  sdft_roi_t roi;
  sdft_fdx_t* twiddles;
};

typedef struct sdft_plan_synthesis sdft_synthesis_t;

struct sdft_plan
{
  sdft_size_t dftsize;
  sdft_analysis_t analysis;
  sdft_synthesis_t synthesis;
};

typedef struct sdft_plan sdft_t;

const sdft_size_t sdft_kernel_size = 2;

sdft_td_t sdft_etc_exchange(sdft_td_t* oldvalue, const sdft_td_t newvalue)
{
  const sdft_td_t value = *oldvalue;
  *oldvalue = newvalue;
  return value;
}

sdft_fd_t sdft_etc_cos(const sdft_fd_t arg)
{
  #if defined(SDFT_FD_FLOAT)
    return cosf(arg);
  #elif defined(SDFT_FD_DOUBLE)
    return cos(arg);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return cosl(arg);
  #endif
}

sdft_fd_t sdft_etc_acos(const sdft_fd_t arg)
{
  #if defined(SDFT_FD_FLOAT)
    return acosf(arg);
  #elif defined(SDFT_FD_DOUBLE)
    return acos(arg);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return acosl(arg);
  #endif
}

sdft_fd_t sdft_etc_real(const sdft_fdx_t z)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return z.r;
  #else
    #if defined(SDFT_FD_FLOAT)
      return crealf(z);
    #elif defined(SDFT_FD_DOUBLE)
      return creal(z);
    #elif defined(SDFT_FD_LONG_DOUBLE)
      return creall(z);
    #endif
  #endif
}

sdft_fd_t sdft_etc_imag(const sdft_fdx_t z)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return z.i;
  #else
    #if defined(SDFT_FD_FLOAT)
      return cimagf(z);
    #elif defined(SDFT_FD_DOUBLE)
      return cimag(z);
    #elif defined(SDFT_FD_LONG_DOUBLE)
      return cimagl(z);
    #endif
  #endif
}

sdft_fdx_t sdft_etc_complex(const sdft_fd_t r, const sdft_fd_t i)
{
  return (sdft_fdx_t){ r, i };
}

sdft_fdx_t sdft_etc_conj(const sdft_fdx_t z)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return sdft_etc_complex(z.r, -z.i);
  #else
    #if defined(SDFT_FD_FLOAT)
      return conjf(z);
    #elif defined(SDFT_FD_DOUBLE)
      return conj(z);
    #elif defined(SDFT_FD_LONG_DOUBLE)
      return conjl(z);
    #endif
  #endif
}

sdft_fdx_t sdft_etc_add(const sdft_fdx_t x, const sdft_fdx_t y)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return sdft_etc_complex(
      x.r + y.r,
      x.i + y.i);
  #elif defined(SDFT_MSVC)
    return sdft_etc_complex(
      sdft_etc_real(x) + sdft_etc_real(y),
      sdft_etc_imag(x) + sdft_etc_imag(y));
  #else
    return x + y;
  #endif
}

sdft_fdx_t sdft_etc_sub(const sdft_fdx_t x, const sdft_fdx_t y)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return sdft_etc_complex(
      x.r - y.r,
      x.i - y.i);
  #elif defined(SDFT_MSVC)
    return sdft_etc_complex(
      sdft_etc_real(x) - sdft_etc_real(y),
      sdft_etc_imag(x) - sdft_etc_imag(y));
  #else
    return x - y;
  #endif
}

sdft_fdx_t sdft_etc_mul(const sdft_fdx_t x, const sdft_fdx_t y)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return sdft_etc_complex(
      x.r * y.r - x.i * y.i,
      x.r * y.i + x.i * y.r);
  #elif defined(SDFT_MSVC)
    #if defined(SDFT_FD_FLOAT)
      return _FCmulcc(x, y);
    #elif defined(SDFT_FD_DOUBLE)
      return _Cmulcc(x, y);
    #elif defined(SDFT_FD_LONG_DOUBLE)
      return _LCmulcc(x, y);
    #endif
  #else
    return x * y;
  #endif
}

sdft_fdx_t sdft_etc_mul_real(const sdft_fd_t x, const sdft_fdx_t y)
{
  #if defined(SDFT_NO_COMPLEX_H)
    return sdft_etc_complex(
      x * y.r,
      x * y.i);
  #elif defined(SDFT_MSVC)
    #if defined(SDFT_FD_FLOAT)
      return _FCmulcr(y, x);
    #elif defined(SDFT_FD_DOUBLE)
      return _Cmulcr(y, x);
    #elif defined(SDFT_FD_LONG_DOUBLE)
      return _LCmulcr(y, x);
    #endif
  #else
    return x * y;
  #endif
}

sdft_fdx_t sdft_etc_polar(const sdft_fd_t r, const sdft_fd_t t)
{
  #if defined(SDFT_FD_FLOAT)
    return sdft_etc_complex(
      r * cosf(t),
      r * sinf(t));
  #elif defined(SDFT_FD_DOUBLE)
    return sdft_etc_complex(
      r * cos(t),
      r * sin(t));
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return sdft_etc_complex(
      r * cosl(t),
      r * sinl(t));
  #endif
}

sdft_fdx_t sdft_etc_convolve(const sdft_fdx_t* values, const sdft_window_t window, const sdft_fd_t weight)
{
  const sdft_size_t l2 = sdft_kernel_size - 2;
  const sdft_size_t l1 = sdft_kernel_size - 1;
  const sdft_size_t m  = sdft_kernel_size;
  const sdft_size_t r1 = sdft_kernel_size + 1;
  const sdft_size_t r2 = sdft_kernel_size + 2;

  switch (window)
  {
    case sdft_window_hann:
    {
      const sdft_fdx_t a = sdft_etc_add(values[m], values[m]);
      const sdft_fdx_t b = sdft_etc_add(values[l1], values[r1]);

      return sdft_etc_mul_real((sdft_fd_t)(0.25) * weight, sdft_etc_sub(a, b));
    }
    case sdft_window_hamming:
    {
      const sdft_fdx_t a = sdft_etc_mul_real((sdft_fd_t)(0.54), values[m]);
      const sdft_fdx_t b = sdft_etc_mul_real((sdft_fd_t)(0.23), sdft_etc_add(values[l1], values[r1]));

      return sdft_etc_mul_real(weight, sdft_etc_sub(a, b));
    }
    case sdft_window_blackman:
    {
      const sdft_fdx_t a = sdft_etc_mul_real((sdft_fd_t)(0.42), values[m]);
      const sdft_fdx_t b = sdft_etc_mul_real((sdft_fd_t)(0.25), sdft_etc_add(values[l1], values[r1]));
      const sdft_fdx_t c = sdft_etc_mul_real((sdft_fd_t)(0.04), sdft_etc_add(values[l2], values[r2]));

      return sdft_etc_mul_real(weight, sdft_etc_add(sdft_etc_sub(a, b), c));
    }
    default:
    {
      return sdft_etc_mul_real(weight, values[m]);
    }
  }
}

/**
 * Allocates a new SDFT plan.
 * @param dftsize Desired number of DFT bins.
 * @param window Analysis window type (boxcar, hann, hamming or blackman).
 * @param latency Synthesis latency factor between 0 and 1.
 *                The default value 1 corresponds to the highest latency and best possible SNR.
 *                A smaller value decreases both latency and SNR, but also increases the workload.
 * @return SDFT plan instance.
 **/
sdft_t* sdft_alloc_custom(const sdft_size_t dftsize, const sdft_window_t window, const sdft_double_t latency)
{
  sdft_t* sdft = (sdft_t*)malloc(sizeof(sdft_t));

  sdft->dftsize = dftsize;

  sdft->analysis.window = window;
  sdft->synthesis.latency = latency;

  sdft->analysis.weight = (sdft_fd_t)(1) / (dftsize * 2);
  sdft->synthesis.weight = (sdft_fd_t)(2);

  sdft->analysis.roi = (sdft_roi_t){ 0, dftsize };
  sdft->synthesis.roi = (sdft_roi_t){ 0, dftsize };

  sdft->analysis.twiddles = (sdft_fdx_t*)calloc(dftsize, sizeof(sdft_fdx_t));
  sdft->synthesis.twiddles = (sdft_fdx_t*)calloc(dftsize, sizeof(sdft_fdx_t));

  sdft->analysis.cursor = 0;
  sdft->analysis.maxcursor = dftsize * 2 - 1;
  sdft->analysis.input = (sdft_td_t*)calloc(dftsize * 2, sizeof(sdft_td_t));

  sdft->analysis.accoutput = (sdft_fdx_t*)calloc(dftsize, sizeof(sdft_fdx_t));
  sdft->analysis.auxoutput = (sdft_fdx_t*)calloc(dftsize + sdft_kernel_size * 2, sizeof(sdft_fdx_t));
  sdft->analysis.fiddles = (sdft_fdx_t*)calloc(dftsize, sizeof(sdft_fdx_t));

  const sdft_fd_t omega = (sdft_fd_t)(-2) * sdft_etc_acos((sdft_fd_t)(-1)) / (dftsize * 2);
  const sdft_fd_t weight = (sdft_fd_t)(+2) / ((sdft_fd_t)(1) - sdft_etc_cos(omega * dftsize * latency));

  for (sdft_size_t i = 0; i < dftsize; ++i)
  {
    sdft->analysis.twiddles[i] = sdft_etc_polar((sdft_fd_t)(1), omega * i);
    sdft->synthesis.twiddles[i] = sdft_etc_polar(weight, omega * i * dftsize * latency);
    sdft->analysis.fiddles[i] = sdft_etc_complex(1, 0);
  }

  return sdft;
}

/**
 * Allocates a new SDFT plan.
 * @param dftsize Desired number of DFT bins.
 * @return SDFT plan instance.
 **/
sdft_t* sdft_alloc(const sdft_size_t dftsize)
{
  return sdft_alloc_custom(dftsize, sdft_window_hann, 1);
}

/**
 * Releases the allocated SDFT plan.
 * @param sdft SDFT plan instance.
 **/
void sdft_free(sdft_t* sdft)
{
  if (sdft == NULL)
  {
    return;
  }

  if (sdft->analysis.twiddles != NULL)
  {
    free(sdft->analysis.twiddles);
    sdft->analysis.twiddles = NULL;
  }

  if (sdft->synthesis.twiddles != NULL)
  {
    free(sdft->synthesis.twiddles);
    sdft->synthesis.twiddles = NULL;
  }

  if (sdft->analysis.input != NULL)
  {
    free(sdft->analysis.input);
    sdft->analysis.input = NULL;
  }

  if (sdft->analysis.accoutput != NULL)
  {
    free(sdft->analysis.accoutput);
    sdft->analysis.accoutput = NULL;
  }

  if (sdft->analysis.auxoutput != NULL)
  {
    free(sdft->analysis.auxoutput);
    sdft->analysis.auxoutput = NULL;
  }

  if (sdft->analysis.fiddles != NULL)
  {
    free(sdft->analysis.fiddles);
    sdft->analysis.fiddles = NULL;
  }

  free(sdft);
  sdft = NULL;
}

/**
 * Resets the specified SDFT plan instance to its initial state.
 * @param sdft SDFT plan instance.
 **/
void sdft_reset(sdft_t* sdft)
{
  sdft->analysis.cursor = 0;

  memset(sdft->analysis.input, 0, (sdft->dftsize * 2) * sizeof(sdft_td_t));
  memset(sdft->analysis.accoutput, 0, (sdft->dftsize) * sizeof(sdft_fdx_t));
  memset(sdft->analysis.auxoutput, 0, (sdft->dftsize + sdft_kernel_size * 2) * sizeof(sdft_fdx_t));

  for (sdft_size_t i = 0; i < sdft->dftsize; ++i)
  {
    sdft->analysis.fiddles[i] = sdft_etc_complex(1, 0);
  }
}

/**
 * Returns the assigned number of DFT bins.
 * @param sdft SDFT plan instance.
 **/
sdft_size_t sdft_size(const sdft_t* sdft)
{
  return (sdft != NULL) ? sdft->dftsize : 0;
}

/**
 * Returns the assigned analysis window type.
 **/
sdft_window_t sdft_window(const sdft_t* sdft)
{
  return (sdft != NULL) ? sdft->analysis.window : sdft_window_boxcar;
}

/**
 * Returns the assigned synthesis latency factor.
 **/
sdft_double_t sdft_latency(const sdft_t* sdft)
{
  return (sdft != NULL) ? sdft->synthesis.latency : 0;
}

/**
 * Estimates the DFT vector for the given sample.
 * @param sdft SDFT plan instance.
 * @param sample Single sample to be analyzed.
 * @param dft Already allocated DFT vector of shape (dftsize).
 **/
void sdft_sdft(sdft_t* sdft, const sdft_td_t sample, sdft_fdx_t* const dft)
{
  const sdft_fd_t delta = sample - sdft_etc_exchange(&sdft->analysis.input[sdft->analysis.cursor], sample);

  if (sdft->analysis.cursor >= sdft->analysis.maxcursor)
  {
    sdft->analysis.cursor = 0;

    for (sdft_size_t i = sdft->analysis.roi.first, j = i + sdft_kernel_size; i < sdft->analysis.roi.second; ++i, ++j)
    {
      sdft->analysis.accoutput[i] = sdft_etc_add(sdft->analysis.accoutput[i], sdft_etc_mul_real(delta, sdft->analysis.fiddles[i]));
      sdft->analysis.fiddles[i]   = sdft_etc_complex(1, 0);
      sdft->analysis.auxoutput[j] = sdft->analysis.accoutput[i];
    }
  }
  else
  {
    sdft->analysis.cursor += 1;

    for (sdft_size_t i = sdft->analysis.roi.first, j = i + sdft_kernel_size; i < sdft->analysis.roi.second; ++i, ++j)
    {
      sdft->analysis.accoutput[i] = sdft_etc_add(sdft->analysis.accoutput[i], sdft_etc_mul_real(delta, sdft->analysis.fiddles[i]));
      sdft->analysis.fiddles[i]   = sdft_etc_mul(sdft->analysis.fiddles[i], sdft->analysis.twiddles[i]);
      sdft->analysis.auxoutput[j] = sdft_etc_mul(sdft->analysis.accoutput[i], sdft_etc_conj(sdft->analysis.fiddles[i]));
    }
  }

  const sdft_size_t auxoffset[] = { sdft_kernel_size, sdft_kernel_size + (sdft->dftsize - 1) };

  for (sdft_size_t i = 1; i <= sdft_kernel_size; ++i)
  {
    sdft->analysis.auxoutput[auxoffset[0] - i] = sdft_etc_conj(sdft->analysis.auxoutput[auxoffset[0] + i]);
    sdft->analysis.auxoutput[auxoffset[1] + i] = sdft_etc_conj(sdft->analysis.auxoutput[auxoffset[1] - i]);
  }

  for (sdft_size_t i = sdft->analysis.roi.first; i < sdft->analysis.roi.second; ++i)
  {
    dft[i] = sdft_etc_convolve(sdft->analysis.auxoutput + i, sdft->analysis.window, sdft->analysis.weight);
  }
}

/**
 * Estimates the DFT matrix for the given sample array.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be analyzed.
 * @param samples Sample array of shape (nsamples).
 * @param dfts Already allocated DFT matrix of shape (nsamples, dftsize).
 **/
void sdft_sdft_n(sdft_t* sdft, const sdft_size_t nsamples, const sdft_td_t* samples, sdft_fdx_t* const dfts)
{
  for (sdft_size_t i = 0; i < nsamples; ++i)
  {
    sdft_sdft(sdft, samples[i], &dfts[i * sdft->dftsize]);
  }
}

/**
 * Estimates the DFT matrix for the given sample array.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be analyzed.
 * @param samples Sample array of shape (nsamples).
 * @param dfts Already allocated array of DFT vectors of shape (nsamples) and (dftsize) correspondingly.
 **/
void sdft_sdft_nd(sdft_t* sdft, const sdft_size_t nsamples, const sdft_td_t* samples, sdft_fdx_t** const dfts)
{
  for (sdft_size_t i = 0; i < nsamples; ++i)
  {
    sdft_sdft(sdft, samples[i], dfts[i]);
  }
}

/**
 * Synthesizes a single sample from the given DFT vector.
 * @param sdft SDFT plan instance.
 * @param dft DFT vector of shape (dftsize).
 **/
sdft_td_t sdft_isdft(sdft_t* sdft, const sdft_fdx_t* dft)
{
  sdft_fd_t sample = (sdft_fd_t)(0);

  if (sdft->synthesis.latency == 1)
  {
    for (sdft_size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_etc_real(dft[i]) * (i % 2 ? -1 : +1);
    }
  }
  else
  {
    for (sdft_size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_etc_real(sdft_etc_mul(dft[i], sdft->synthesis.twiddles[i]));
    }
  }

  sample *= sdft->synthesis.weight;

  return (sdft_td_t)(sample);
}

/**
 * Synthesizes the sample array from the given DFT matrix.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be synthesized.
 * @param dfts DFT matrix of shape (nsamples, dftsize).
 * @param samples Already allocated sample array of shape (nsamples).
 **/
void sdft_isdft_n(sdft_t* sdft, const sdft_size_t nsamples, const sdft_fdx_t* dfts, sdft_td_t* const samples)
{
  for (sdft_size_t i = 0; i < nsamples; ++i)
  {
    samples[i] = sdft_isdft(sdft, &dfts[i * sdft->dftsize]);
  }
}

/**
 * Synthesizes the sample array from the given DFT matrix.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be synthesized.
 * @param dfts Array of DFT vectors of shape (nsamples) and (dftsize) correspondingly.
 * @param samples Already allocated sample array of shape (nsamples).
 **/
void sdft_isdft_nd(sdft_t* sdft, const sdft_size_t nsamples, const sdft_fdx_t** dfts, sdft_td_t* const samples)
{
  for (sdft_size_t i = 0; i < nsamples; ++i)
  {
    samples[i] = sdft_isdft(sdft, dfts[i]);
  }
}

#if defined(__cplusplus)
}
#endif
