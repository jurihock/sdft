#ifndef SDFT_H
#define SDFT_H

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#if defined(SDFT_TD_FLOAT)
  typedef float SDFT_TD_TYPE;
#elif defined(SDFT_TD_DOUBLE)
  typedef double SDFT_TD_TYPE;
#elif defined(SDFT_TD_LONG_DOUBLE)
  typedef long double SDFT_TD_TYPE;
#else
  #define SDFT_TD_FLOAT
  typedef float SDFT_TD_TYPE;
#endif

#if defined(SDFT_FD_FLOAT)
  typedef float SDFT_FD_TYPE;
  typedef float complex SDFT_FDX_TYPE;
#elif defined(SDFT_FD_DOUBLE)
  typedef double SDFT_FD_TYPE;
  typedef double complex SDFT_FDX_TYPE;
#elif defined(SDFT_FD_LONG_DOUBLE)
  typedef long double SDFT_FD_TYPE;
  typedef long double complex SDFT_FDX_TYPE;
#else
  #define SDFT_FD_DOUBLE
  typedef double SDFT_FD_TYPE;
  typedef double complex SDFT_FDX_TYPE;
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct sdft_plan_roi
{
  size_t first;
  size_t second;
};

typedef struct sdft_plan_roi SDFT_ROI;

struct sdft_plan_analysis
{
  SDFT_ROI roi;
  SDFT_FDX_TYPE* twiddles;

  size_t cursor;
  size_t maxcursor;
  SDFT_TD_TYPE* input;
  SDFT_FDX_TYPE* accoutput;
  SDFT_FDX_TYPE* auxoutput;
  SDFT_FDX_TYPE* fiddles;
};

typedef struct sdft_plan_analysis SDFT_ANALYSIS;

struct sdft_plan_synthesis
{
  SDFT_ROI roi;
  SDFT_FDX_TYPE* twiddles;
};

typedef struct sdft_plan_synthesis SDFT_SYNTHESIS;

struct sdft_plan
{
  size_t dftsize;
  double latency;

  SDFT_ANALYSIS analysis;
  SDFT_SYNTHESIS synthesis;
};

typedef struct sdft_plan SDFT;

#if defined(SDFT_FD_FLOAT)
  SDFT_FD_TYPE sdft_acos(const SDFT_FD_TYPE a) { return acosf(a); }
#elif defined(SDFT_FD_DOUBLE)
  SDFT_FD_TYPE sdft_acos(const SDFT_FD_TYPE a) { return acos(a); }
#elif defined(SDFT_FD_LONG_DOUBLE)
  SDFT_FD_TYPE sdft_acos(const SDFT_FD_TYPE a) { return acosl(a); }
#endif

#if defined(SDFT_FD_FLOAT)
  SDFT_FD_TYPE sdft_cos(const SDFT_FD_TYPE a) { return cosf(a); }
#elif defined(SDFT_FD_DOUBLE)
  SDFT_FD_TYPE sdft_cos(const SDFT_FD_TYPE a) { return cos(a); }
#elif defined(SDFT_FD_LONG_DOUBLE)
  SDFT_FD_TYPE sdft_cos(const SDFT_FD_TYPE a) { return cosl(a); }
#endif

#if defined(SDFT_FD_FLOAT)
  SDFT_FD_TYPE sdft_real(const SDFT_FDX_TYPE z) { return crealf(z); }
#elif defined(SDFT_FD_DOUBLE)
  SDFT_FD_TYPE sdft_real(const SDFT_FDX_TYPE z) { return creal(z); }
#elif defined(SDFT_FD_LONG_DOUBLE)
  SDFT_FD_TYPE sdft_real(const SDFT_FDX_TYPE z) { return creall(z); }
#endif

#if defined(SDFT_FD_FLOAT)
  SDFT_FDX_TYPE sdft_polar(const SDFT_FD_TYPE r, const SDFT_FD_TYPE a) { return r * cexpf(I * a); }
#elif defined(SDFT_FD_DOUBLE)
  SDFT_FDX_TYPE sdft_polar(const SDFT_FD_TYPE r, const SDFT_FD_TYPE a) { return r * cexp(I * a); }
#elif defined(SDFT_FD_LONG_DOUBLE)
  SDFT_FDX_TYPE sdft_polar(const SDFT_FD_TYPE r, const SDFT_FD_TYPE a) { return r * cexpl(I * a); }
#endif

#if defined(SDFT_FD_FLOAT)
  SDFT_FDX_TYPE sdft_conj(const SDFT_FDX_TYPE z) { return conjf(z); }
#elif defined(SDFT_FD_DOUBLE)
  SDFT_FDX_TYPE sdft_conj(const SDFT_FDX_TYPE z) { return conj(z); }
#elif defined(SDFT_FD_LONG_DOUBLE)
  SDFT_FDX_TYPE sdft_conj(const SDFT_FDX_TYPE z) { return conjl(z); }
#endif

SDFT_TD_TYPE sdft_exchange(SDFT_TD_TYPE* old_value, const SDFT_TD_TYPE new_value)
{
  SDFT_TD_TYPE value = *old_value;
  *old_value = new_value;
  return value;
}

SDFT_FDX_TYPE sdft_window(const SDFT_FDX_TYPE left,
                          const SDFT_FDX_TYPE middle,
                          const SDFT_FDX_TYPE right,
                          const SDFT_FD_TYPE  weight)
{
  return (SDFT_FD_TYPE)(0.25) * ((middle + middle) - (left + right)) * weight;
}

SDFT* sdft_alloc_custom(const size_t dftsize, const double latency)
{
  SDFT* sdft = malloc(sizeof(SDFT));

  assert(sdft != NULL);

  sdft->dftsize = dftsize;
  sdft->latency = latency;

  sdft->analysis.roi = (SDFT_ROI){ 0, dftsize };
  sdft->synthesis.roi = (SDFT_ROI){ 0, dftsize };

  sdft->analysis.twiddles = calloc(dftsize, sizeof(SDFT_FDX_TYPE));
  sdft->synthesis.twiddles = calloc(dftsize, sizeof(SDFT_FDX_TYPE));

  assert(sdft->analysis.twiddles != NULL);
  assert(sdft->synthesis.twiddles != NULL);

  sdft->analysis.cursor = 0;
  sdft->analysis.maxcursor = dftsize * 2 - 1;
  sdft->analysis.input = calloc(dftsize * 2, sizeof(SDFT_TD_TYPE));
  sdft->analysis.accoutput = calloc(dftsize, sizeof(SDFT_FDX_TYPE));
  sdft->analysis.auxoutput = calloc(dftsize + 2, sizeof(SDFT_FDX_TYPE));
  sdft->analysis.fiddles = calloc(dftsize, sizeof(SDFT_FDX_TYPE));

  assert(sdft->analysis.input != NULL);
  assert(sdft->analysis.accoutput != NULL);
  assert(sdft->analysis.auxoutput != NULL);
  assert(sdft->analysis.fiddles != NULL);

  const SDFT_FD_TYPE pi = (SDFT_FD_TYPE)(-2) * sdft_acos((SDFT_FD_TYPE)(-1)) / (dftsize * 2);
  const SDFT_FD_TYPE weight = (SDFT_FD_TYPE)(2) / ((SDFT_FD_TYPE)(1) - sdft_cos(pi * dftsize * latency));

  for (size_t i = 0; i < dftsize; ++i)
  {
    sdft->analysis.twiddles[i] = sdft_polar((SDFT_FD_TYPE)(1), pi * i);
    sdft->synthesis.twiddles[i] = sdft_polar(weight, pi * i * dftsize * latency);

    sdft->analysis.fiddles[i] = (SDFT_FD_TYPE)(1);
  }

  return sdft;
}

SDFT* sdft_alloc(const size_t dftsize)
{
  return sdft_alloc_custom(dftsize, 1);
}

void sdft_free(SDFT* sdft)
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

size_t sdft_size(const SDFT* sdft)
{
  return (sdft != NULL) ? sdft->dftsize : 0;
}

void sdft_sdft(SDFT* sdft, const SDFT_TD_TYPE sample, SDFT_FDX_TYPE* const dft)
{
  // assert(dft.size() == sdft->dftsize);

  // actually the weight denominator needs to be dftsize*2 to get proper magnitude scaling,
  // but then requires a correction by factor 2 in synthesis and is therefore omitted

  const SDFT_FD_TYPE weight = (SDFT_FD_TYPE)(1) / sdft->dftsize;

  const SDFT_FD_TYPE delta = sample - sdft_exchange(&sdft->analysis.input[sdft->analysis.cursor], sample);

  for (size_t i = sdft->analysis.roi.first, j = i + 1; i < sdft->analysis.roi.second; ++i, ++j)
  {
    const SDFT_FDX_TYPE oldfiddle = sdft->analysis.fiddles[i];
    const SDFT_FDX_TYPE newfiddle = oldfiddle * sdft->analysis.twiddles[i];

    sdft->analysis.fiddles[i] = newfiddle;

    sdft->analysis.accoutput[i] += delta * oldfiddle;
    sdft->analysis.auxoutput[j] = sdft->analysis.accoutput[i] * sdft_conj(newfiddle);
  }

  // theoretically the DFT periodicity needs to be preserved for proper windowing,
  // but the both outer bins seem to be noisy for an unclear reason
  // and will be suppressed anyway after windowing

  // analysis.auxoutput[0] = analysis.auxoutput[sdft->dftsize];
  // analysis.auxoutput[sdft->dftsize + 1] = analysis.auxoutput[1];

  for (size_t i = sdft->analysis.roi.first, j = i + 1; i < sdft->analysis.roi.second; ++i, ++j)
  {
    dft[i] = sdft_window(sdft->analysis.auxoutput[j - 1],
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

SDFT_TD_TYPE sdft_isdft(SDFT* sdft, const SDFT_FDX_TYPE* dft)
{
  // assert(dft.size() == dftsize);

  SDFT_FD_TYPE sample = (SDFT_FD_TYPE)(0);

  if (sdft->latency == 1)
  {
    for (size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_real(dft[i]) * (i % 2 ? -1 : +1);
    }
  }
  else
  {
    for (size_t i = sdft->synthesis.roi.first; i < sdft->synthesis.roi.second; ++i)
    {
      sample += sdft_real(dft[i] * sdft->synthesis.twiddles[i]);
    }
  }

  return (SDFT_TD_TYPE)(sample);
}

#ifdef __cplusplus
}
#endif

#endif
