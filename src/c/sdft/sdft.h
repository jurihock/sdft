/**
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
 * #define SDFT_FD_DOUBLE       // default, recommended
 * #define SDFT_FD_LONG_DOUBLE
 *
 * The specified data type appears as typedef
 * sdft_td_t, sdft_fd_t and sdft_fdx_t for complex numbers.
 **/

/**
 * List of common functions:
 *
 * sdft_alloc
 * sdft_free
 *
 * sdft_size
 *
 * sdft_sdft
 * sdft_sdft_n
 * sdft_sdft_nd
 *
 * sdft_isdft
 * sdft_isdft_n
 * sdft_isdft_nd
 **/

#pragma once

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(SDFT_TD_FLOAT)
  typedef float sdft_td_t;
#elif defined(SDFT_TD_DOUBLE)
  typedef double sdft_td_t;
#elif defined(SDFT_TD_LONG_DOUBLE)
  typedef long double sdft_td_t;
#else
  #define SDFT_TD_FLOAT
  typedef float sdft_td_t;
#endif

#if defined(SDFT_FD_FLOAT)
  typedef float sdft_fd_t;
  typedef float complex sdft_fdx_t;
#elif defined(SDFT_FD_DOUBLE)
  typedef double sdft_fd_t;
  typedef double complex sdft_fdx_t;
#elif defined(SDFT_FD_LONG_DOUBLE)
  typedef long double sdft_fd_t;
  typedef long double complex sdft_fdx_t;
#else
  #define SDFT_FD_DOUBLE
  typedef double sdft_fd_t;
  typedef double complex sdft_fdx_t;
#endif

struct sdft_plan_roi
{
  size_t first;
  size_t second;
};

typedef struct sdft_plan_roi sdft_roi_t;

struct sdft_plan_analysis
{
  sdft_roi_t roi;
  sdft_fdx_t* twiddles;

  size_t cursor;
  size_t maxcursor;
  sdft_td_t* input;

  sdft_fdx_t* accoutput;
  sdft_fdx_t* auxoutput;
  sdft_fdx_t* fiddles;
};

typedef struct sdft_plan_analysis sdft_analysis_t;

struct sdft_plan_synthesis
{
  sdft_roi_t roi;
  sdft_fdx_t* twiddles;
};

typedef struct sdft_plan_synthesis sdft_synthesis_t;

struct sdft_plan
{
  size_t dftsize;
  double latency;

  sdft_analysis_t analysis;
  sdft_synthesis_t synthesis;
};

typedef struct sdft_plan sdft_t;

sdft_fd_t sdft_etc_acos(const sdft_fd_t a)
{
  #if defined(SDFT_FD_FLOAT)
    return acosf(a);
  #elif defined(SDFT_FD_DOUBLE)
    return acos(a);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return acosl(a);
  #endif
}

sdft_fd_t sdft_etc_cos(const sdft_fd_t a)
{
  #if defined(SDFT_FD_FLOAT)
    return cosf(a);
  #elif defined(SDFT_FD_DOUBLE)
    return cos(a);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return cosl(a);
  #endif
}

sdft_fd_t sdft_etc_real(const sdft_fdx_t z)
{
  #if defined(SDFT_FD_FLOAT)
    return crealf(z);
  #elif defined(SDFT_FD_DOUBLE)
    return creal(z);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return creall(z);
  #endif
}

sdft_fdx_t sdft_etc_polar(const sdft_fd_t r, const sdft_fd_t a)
{
  #if defined(SDFT_FD_FLOAT)
    return r * cexpf(I * a);
  #elif defined(SDFT_FD_DOUBLE)
    return r * cexp(I * a);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return r * cexpl(I * a);
  #endif
}

sdft_fdx_t sdft_etc_conj(const sdft_fdx_t z)
{
  #if defined(SDFT_FD_FLOAT)
    return conjf(z);
  #elif defined(SDFT_FD_DOUBLE)
    return conj(z);
  #elif defined(SDFT_FD_LONG_DOUBLE)
    return conjl(z);
  #endif
}

sdft_td_t sdft_etc_exchange(sdft_td_t* old_value, const sdft_td_t new_value)
{
  sdft_td_t value = *old_value;
  *old_value = new_value;
  return value;
}

sdft_fdx_t sdft_etc_window(const sdft_fdx_t left, const sdft_fdx_t middle, const sdft_fdx_t right, const sdft_fd_t  weight)
{
  return (sdft_fd_t)(0.25) * ((middle + middle) - (left + right)) * weight;
}

/**
 * Allocates a new SDFT plan.
 * @param dftsize Desired number of DFT bins.
 * @param latency Synthesis latency factor between 0 and 1.
 *                The default value 1 corresponds to the highest latency and best possible SNR.
 *                A smaller value decreases both latency and SNR, but also increases the workload.
 * @return SDFT plan instance.
 **/
sdft_t* sdft_alloc_custom(const size_t dftsize, const double latency)
{
  sdft_t* sdft = malloc(sizeof(sdft_t));

  assert(sdft != NULL);

  sdft->dftsize = dftsize;
  sdft->latency = latency;

  sdft->analysis.roi = (sdft_roi_t){ 0, dftsize };
  sdft->synthesis.roi = (sdft_roi_t){ 0, dftsize };

  sdft->analysis.twiddles = calloc(dftsize, sizeof(sdft_fdx_t));
  sdft->synthesis.twiddles = calloc(dftsize, sizeof(sdft_fdx_t));

  assert(sdft->analysis.twiddles != NULL);
  assert(sdft->synthesis.twiddles != NULL);

  sdft->analysis.cursor = 0;
  sdft->analysis.maxcursor = dftsize * 2 - 1;
  sdft->analysis.input = calloc(dftsize * 2, sizeof(sdft_td_t));

  sdft->analysis.accoutput = calloc(dftsize, sizeof(sdft_fdx_t));
  sdft->analysis.auxoutput = calloc(dftsize + 2, sizeof(sdft_fdx_t));
  sdft->analysis.fiddles = calloc(dftsize, sizeof(sdft_fdx_t));

  assert(sdft->analysis.input != NULL);
  assert(sdft->analysis.accoutput != NULL);
  assert(sdft->analysis.auxoutput != NULL);
  assert(sdft->analysis.fiddles != NULL);

  const sdft_fd_t pi = (sdft_fd_t)(-2) * sdft_etc_acos((sdft_fd_t)(-1)) / (dftsize * 2);
  const sdft_fd_t weight = (sdft_fd_t)(2) / ((sdft_fd_t)(1) - sdft_etc_cos(pi * dftsize * latency));

  for (size_t i = 0; i < dftsize; ++i)
  {
    sdft->analysis.twiddles[i] = sdft_etc_polar((sdft_fd_t)(1), pi * i);
    sdft->synthesis.twiddles[i] = sdft_etc_polar(weight, pi * i * dftsize * latency);

    sdft->analysis.fiddles[i] = (sdft_fd_t)(1);
  }

  return sdft;
}

/**
 * Allocates a new SDFT plan.
 * @param dftsize Desired number of DFT bins.
 * @return SDFT plan instance.
 **/
sdft_t* sdft_alloc(const size_t dftsize)
{
  return sdft_alloc_custom(dftsize, 1);
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
 * Returns the assigned number of DFT bins.
 * @param sdft SDFT plan instance.
 **/
size_t sdft_size(const sdft_t* sdft)
{
  return (sdft != NULL) ? sdft->dftsize : 0;
}

/**
 * Estimates the DFT vector for the given sample.
 * @param sdft SDFT plan instance.
 * @param sample Single sample to be analyzed.
 * @param dft Already allocated DFT vector of shape (dftsize).
 **/
void sdft_sdft(sdft_t* sdft, const sdft_td_t sample, sdft_fdx_t* const dft)
{
  // assert(dft.size() == sdft->dftsize);

  // actually the weight denominator needs to be dftsize*2 to get proper magnitude scaling,
  // but then requires a correction by factor 2 in synthesis and is therefore omitted

  const sdft_fd_t weight = (sdft_fd_t)(1) / sdft->dftsize;

  const sdft_fd_t delta = sample - sdft_etc_exchange(&sdft->analysis.input[sdft->analysis.cursor], sample);

  for (size_t i = sdft->analysis.roi.first, j = i + 1; i < sdft->analysis.roi.second; ++i, ++j)
  {
    const sdft_fdx_t oldfiddle = sdft->analysis.fiddles[i];
    const sdft_fdx_t newfiddle = oldfiddle * sdft->analysis.twiddles[i];

    sdft->analysis.fiddles[i] = newfiddle;

    sdft->analysis.accoutput[i] += delta * oldfiddle;
    sdft->analysis.auxoutput[j] = sdft->analysis.accoutput[i] * sdft_etc_conj(newfiddle);
  }

  // theoretically the DFT periodicity needs to be preserved for proper windowing,
  // but the both outer bins seem to be noisy for an unclear reason
  // and will be suppressed anyway after windowing

  // analysis.auxoutput[0] = analysis.auxoutput[sdft->dftsize];
  // analysis.auxoutput[sdft->dftsize + 1] = analysis.auxoutput[1];

  for (size_t i = sdft->analysis.roi.first, j = i + 1; i < sdft->analysis.roi.second; ++i, ++j)
  {
    dft[i] = sdft_etc_window(sdft->analysis.auxoutput[j - 1],
                             sdft->analysis.auxoutput[j],
                             sdft->analysis.auxoutput[j + 1],
                             weight);
  }

  // finally suppress outer DFT bins as announced in the comment above

  dft[0] = dft[sdft->dftsize - 1] = 0;

  if (++sdft->analysis.cursor > sdft->analysis.maxcursor)
  {
    sdft->analysis.cursor = 0;
  }
}

/**
 * Estimates the DFT matrix for the given sample array.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be analyzed.
 * @param samples Sample array of shape (nsamples).
 * @param dfts Already allocated DFT matrix of shape (nsamples, dftsize).
 **/
void sdft_sdft_n(sdft_t* sdft, const size_t nsamples, const sdft_td_t* samples, sdft_fdx_t* const dfts)
{
  for (size_t i = 0; i < nsamples; ++i)
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
void sdft_sdft_nd(sdft_t* sdft, const size_t nsamples, const sdft_td_t* samples, sdft_fdx_t** const dfts)
{
  for (size_t i = 0; i < nsamples; ++i)
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
  // assert(dft.size() == dftsize);

  sdft_fd_t sample = (sdft_fd_t)(0);

  if (sdft->latency == 1)
  {
    for (size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_etc_real(dft[i]) * (i % 2 ? -1 : +1);
    }
  }
  else
  {
    for (size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_etc_real(dft[i] * sdft->synthesis.twiddles[i]);
    }
  }

  return (sdft_td_t)(sample);
}

/**
 * Synthesizes the sample array from the given DFT matrix.
 * @param sdft SDFT plan instance.
 * @param nsamples Number of samples to be synthesized.
 * @param dfts DFT matrix of shape (nsamples, dftsize).
 * @param samples Already allocated sample array of shape (nsamples).
 **/
void sdft_isdft_n(sdft_t* sdft, const size_t nsamples, const sdft_fdx_t* dfts, sdft_td_t* const samples)
{
  // assert(samples.size() == dfts.size());

  for (size_t i = 0; i < nsamples; ++i)
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
void sdft_isdft_nd(sdft_t* sdft, const size_t nsamples, const sdft_fdx_t** dfts, sdft_td_t* const samples)
{
  // assert(samples.size() == dfts.size());

  for (size_t i = 0; i < nsamples; ++i)
  {
    samples[i] = sdft_isdft(sdft, dfts[i]);
  }
}

#if defined(__cplusplus)
}
#endif
