/**
 * Copyright (c) 2023 Juergen Hock
 *
 * SPDX-License-Identifier: MIT
 *
 * Sliding DFT algorithm implementation according to [1] combined with [2].
 *
 * [1] Richard Lyons
 *     A Fast Guaranteed-Stable Sliding DFT Algorithm
 *     DSPRelated.com (2023)
 *     https://www.dsprelated.com/showarticle/1533.php
 *
 * [2] Russell Bradford and Richard Dobson and John ffitch
 *     Sliding is Smoother than Jumping
 *     International Computer Music Conference (2005)
 *     http://hdl.handle.net/2027/spo.bbp2372.2005.086
 *
 * Source: https://github.com/jurihock/sdft
 **/

#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>
#include <vector>

/**
 * Sliding Discrete Fourier Transform (SDFT).
 * @tparam T Time domain data type, which can be float (default), double or long double.
 * @tparam F Frequency domain data type, which can be float, double (default and recommended) or long double.
 **/
template <typename T = float, typename F = double>
class SDFT
{

public:

  /**
   * Supported SDFT analysis window types.
   **/
  enum class Window
  {
    Boxcar,
    Hann,
    Hamming,
    Blackman
  };

  /**
   * Creates a new SDFT plan.
   * @param dftsize Desired number of DFT bins.
   * @param window Analysis window type (boxcar, hann, hamming or blackman).
   * @param latency Synthesis latency factor between 0 and 1.
   *                The default value 1 corresponds to the highest latency and best possible SNR.
   *                A smaller value decreases both latency and SNR, but also increases the workload.
   **/
  SDFT(const size_t dftsize, const SDFT::Window window = SDFT::Window::Hann, const double latency = 1) :
    dftsize(dftsize)
  {
    analysis.window = window;
    synthesis.latency = latency;

    analysis.weight = F(1) / dftsize;
    synthesis.weight = F(2);

    analysis.roi = { 0, dftsize };
    synthesis.roi = { 0, dftsize };

    analysis.twiddles.resize(dftsize);
    analysis.fiddles.resize(dftsize);
    synthesis.twiddles.resize(dftsize);

    analysis.cursor = 0;
    analysis.inputs.resize(dftsize);
    analysis.resonator.resize(dftsize);
    analysis.outputs.resize(dftsize + kernelsize * 2);

    const F pi = std::acos(F(-1));
    const F omega = F(2) * pi / (dftsize * 2);
    const F weight = F(2) / (F(1) - std::cos(omega * dftsize * latency));

    for (size_t i = 0; i < dftsize; ++i)
    {
      analysis.twiddles[i] = std::polar(F(1), i * omega);
      analysis.fiddles[i] = F(2) * analysis.twiddles[i].real();
      synthesis.twiddles[i] = std::polar(weight, i * omega * dftsize * latency);
    }
  }

  /**
   * Resets this SDFT plan instance to its initial state.
   **/
  void reset()
  {
    analysis.cursor = 0;
    std::fill(analysis.inputs.begin(), analysis.inputs.end(), T(0));
    std::fill(analysis.resonator.begin(), analysis.resonator.end(), std::array<F, 2>());
    std::fill(analysis.outputs.begin(), analysis.outputs.end(), std::complex<F>(F(0)));
  }

  /**
   * Returns the assigned number of DFT bins.
   **/
  size_t size() const
  {
    return dftsize;
  }

  /**
   * Returns the assigned analysis window type.
   **/
  SDFT::Window window() const
  {
    return analysis.window;
  }

  /**
   * Returns the assigned synthesis latency factor.
   **/
  double latency() const
  {
    return synthesis.latency;
  }

  /**
   * Estimates the DFT vector for the given sample.
   * @param sample Single sample to be analyzed.
   * @param dft Already allocated DFT vector of shape (dftsize).
   **/
  void sdft(const T sample, std::complex<F>* const dft)
  {
    const T newsample = sample;
    const T oldsample = exchange(analysis.inputs[analysis.cursor], newsample);

    analysis.cursor = (analysis.cursor + 1) % analysis.inputs.size();

    const T comb[] = { newsample - oldsample, newsample + oldsample };

    for (size_t i = analysis.roi.first, j = i + kernelsize; i < analysis.roi.second; ++i, ++j)
    {
      const F state = comb[i % 2] + analysis.resonator[i][0] * analysis.fiddles[i] - analysis.resonator[i][1];

      analysis.outputs[j] = state * analysis.twiddles[i] - analysis.resonator[i][0];

      analysis.resonator[i][1] = analysis.resonator[i][0];
      analysis.resonator[i][0] = state;
    }

    const size_t offsets[] = { kernelsize, kernelsize + (dftsize - 1) };

    for (size_t i = 1; i <= kernelsize; ++i)
    {
      analysis.outputs[offsets[0] - i] = std::conj(analysis.outputs[offsets[0] + i]);
      analysis.outputs[offsets[1] + i] = std::conj(analysis.outputs[offsets[1] - i]);
    }

    for (size_t i = analysis.roi.first; i < analysis.roi.second; ++i)
    {
      dft[i] = convolve(analysis.outputs.data() + i, analysis.window, analysis.weight);
    }
  }

  /**
   * Estimates the DFT matrix for the given sample array.
   * @param nsamples Number of samples to be analyzed.
   * @param samples Sample array of shape (nsamples).
   * @param dfts Already allocated DFT matrix of shape (nsamples, dftsize).
   **/
  void sdft(const size_t nsamples, const T* samples, std::complex<F>* const dfts)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      sdft(samples[i], &dfts[i * dftsize]);
    }
  }

  /**
   * Estimates the DFT matrix for the given sample array.
   * @param nsamples Number of samples to be analyzed.
   * @param samples Sample array of shape (nsamples).
   * @param dfts Already allocated array of DFT vectors of shape (nsamples) and (dftsize) correspondingly.
   **/
  void sdft(const size_t nsamples, const T* samples, std::complex<F>** const dfts)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      sdft(samples[i], dfts[i]);
    }
  }

  /**
   * Synthesizes a single sample from the given DFT vector.
   * @param dft DFT vector of shape (dftsize).
   **/
  T isdft(const std::complex<F>* dft)
  {
    F sample = F(0);

    if (synthesis.latency == 1)
    {
      for (size_t i = synthesis.roi.first; i < synthesis.roi.second; ++i)
      {
        sample += dft[i].real() * (i % 2 ? -1 : +1);
      }
    }
    else
    {
      for (size_t i = synthesis.roi.first; i < synthesis.roi.second; ++i)
      {
        sample += (dft[i] * synthesis.twiddles[i]).real();
      }
    }

    sample *= synthesis.weight;

    return static_cast<T>(sample);
  }

  /**
   * Synthesizes the sample array from the given DFT matrix.
   * @param nsamples Number of samples to be synthesized.
   * @param dfts DFT matrix of shape (nsamples, dftsize).
   * @param samples Already allocated sample array of shape (nsamples).
   **/
  void isdft(const size_t nsamples, const std::complex<F>* dfts, T* const samples)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      samples[i] = isdft(&dfts[i * dftsize]);
    }
  }

  /**
   * Synthesizes the sample array from the given DFT matrix.
   * @param nsamples Number of samples to be synthesized.
   * @param dfts Array of DFT vectors of shape (nsamples) and (dftsize) correspondingly.
   * @param samples Already allocated sample array of shape (nsamples).
   **/
  void isdft(const size_t nsamples, const std::complex<F>** dfts, T* const samples)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      samples[i] = isdft(dfts[i]);
    }
  }

private:

  static const size_t kernelsize = 2;
  const size_t dftsize;

  struct
  {
    SDFT::Window window;

    F weight;
    std::pair<size_t, size_t> roi;
    std::vector<std::complex<F>> twiddles;
    std::vector<F> fiddles;

    size_t cursor;
    std::vector<T> inputs;
    std::vector<std::array<F, 2>> resonator;
    std::vector<std::complex<F>> outputs;
  }
  analysis;

  struct
  {
    double latency;

    F weight;
    std::pair<size_t, size_t> roi;
    std::vector<std::complex<F>> twiddles;
  }
  synthesis;

  inline static T exchange(T& oldvalue, const T newvalue)
  {
    const T value = oldvalue;
    oldvalue = newvalue;
    return value;
  }

  inline static std::complex<F> convolve(const std::complex<F>* values, const SDFT::Window window, const F weight)
  {
    const auto stride = 2;

    const auto ll = kernelsize - 2 * stride;
    const auto l  = kernelsize - 1 * stride;
    const auto m  = kernelsize;
    const auto r  = kernelsize + 1 * stride;
    const auto rr = kernelsize + 2 * stride;

    switch (window)
    {
      case SDFT::Window::Hann:
      {
        const std::complex<F> a = values[m] + values[m];
        const std::complex<F> b = values[l] + values[r];

        return F(0.25) * weight * (a - b);
      }
      case SDFT::Window::Hamming:
      {
        const std::complex<F> a = F(0.54) * values[m];
        const std::complex<F> b = F(0.23) * (values[l] + values[r]);

        return weight * (a - b);
      }
      case SDFT::Window::Blackman:
      {
        const std::complex<F> a = F(0.42) * values[m];
        const std::complex<F> b = F(0.25) * (values[l] + values[r]);
        const std::complex<F> c = F(0.04) * (values[ll] + values[rr]);

        return weight * (a - b + c);
      }
      default:
      {
        return weight * values[m];
      }
    }
  }

};
