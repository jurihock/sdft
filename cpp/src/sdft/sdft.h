/**
 * Copyright (c) 2025 Juergen Hock
 *
 * SPDX-License-Identifier: MIT
 *
 * Sliding DFT implementation according to [1] combined with [2].
 *
 * [1] Rick Lyons
 *     An Efficient Full-Band Sliding DFT Spectrum Analyzer
 *     DSPRelated.com (2011)
 *     https://www.dsprelated.com/showarticle/1396.php
 *
 * [2] Russell Bradford and Richard Dobson and John ffitch
 *     Sliding is Smoother than Jumping
 *     International Computer Music Conference (2005)
 *     http://hdl.handle.net/2027/spo.bbp2372.2005.086
 *
 * Source: https://github.com/jurihock/sdft
 **/

#pragma once

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>
#include <vector>

namespace sdft
{
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
   * Sliding Discrete Fourier Transform (SDFT).
   * @tparam T Time domain data type, which can be float (default), double or long double.
   * @tparam F Frequency domain data type, which can be float, double (default and recommended) or long double.
   **/
  template <typename T = float, typename F = double>
  class SDFT
  {

  public:

    /**
     * Creates a new SDFT plan.
     * @param samplerate Sample rate in hertz.
     * @param dftsize Desired number of DFT bins.
     * @param window Analysis window type (boxcar, hann, hamming or blackman).
     * @param latency Synthesis latency factor between 0 and 1.
     *                The default value 1 corresponds to the highest latency and best possible SNR.
     *                A smaller value decreases both latency and SNR, but also increases the workload.
     **/
    SDFT(const double samplerate, const size_t dftsize, const Window window = Window::Hann, const double latency = 1)
    {
      config.samplerate = samplerate;
      config.dftsize = dftsize;

      analysis.window = window;
      synthesis.latency = latency;

      analysis.weight = F(1) / (dftsize * 2);
      synthesis.weight = F(2);

      analysis.buffer.resize(dftsize + kernelsize * 2);

      analysis.twiddles.resize(dftsize);
      synthesis.twiddles.resize(dftsize);

      size_t fullsize;

      if (dftsize % 2) // odd dftsize => even fullsize
      {
        fullsize = dftsize * 2 - 2;
        assert(dftsize == (fullsize / 2) + 1);
      }
      else // even dftsize => odd fullsize
      {
        fullsize = dftsize * 2 - 1;
        assert(dftsize == (fullsize + 1) / 2);
      }

      if (fullsize % 2)
      {
        analysis.feedback.first = F(0.5);
        analysis.feedback.last = F(1.0);
        analysis.feedback.all = F(1) / dftsize;
      }
      else
      {
        analysis.feedback.first = F(0.5);
        analysis.feedback.last = F(0.5);
        analysis.feedback.all = F(1) / dftsize;
      }

      const F omega = F(2) * std::acos(F(-1)) / fullsize;
      const F weight = F(2) / (F(1) - std::cos(omega * dftsize * latency));

      for (size_t i = 0; i < dftsize; ++i)
      {
        analysis.twiddles[i] = std::polar(F(1), omega * i);
        synthesis.twiddles[i] = std::polar(weight, omega * i * dftsize * latency); // TODO: -omega
      }

      config.frequencies.resize(dftsize);

      for (size_t i = 0; i < dftsize; ++i)
      {
        config.frequencies[i] = samplerate * i / fullsize;
      }
    }

    /**
     * Resets this SDFT plan instance to its initial state.
     **/
    void reset()
    {
      std::fill(analysis.buffer.begin(), analysis.buffer.end(), 0);
    }

    /**
     * Returns the assigned number of DFT bins.
     **/
    size_t size() const
    {
      return config.dftsize;
    }

    /**
     * Returns the assigned sample rate in hertz.
     **/
    double samplerate() const
    {
      return config.samplerate;
    }

    /**
     * Returns the assigned analysis window type.
     **/
    Window window() const
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
     * Returns the vector of bin frequencies in hertz.
     **/
    const std::vector<double>& frequencies() const
    {
      return config.frequencies;
    }

    /**
     * Estimates the DFT vector for the given sample.
     * @param sample Single sample to be analyzed.
     * @param dft Already allocated DFT vector of shape (dftsize).
     **/
    void sdft(const T sample, std::complex<F>* const dft)
    {
      F feedback = std::real(analysis.buffer[kernelsize]) * analysis.feedback.first
                + std::real(analysis.buffer[kernelsize + (config.dftsize - 1)]) * analysis.feedback.last;

      for (size_t i = 1, j = i + kernelsize; i < (config.dftsize - 1); ++i, ++j)
      {
        feedback += std::real(analysis.buffer[j]);
      }

      const F delta = sample - feedback * analysis.feedback.all;

      for (size_t i = 0, j = i + kernelsize; i < config.dftsize; ++i, ++j)
      {
        analysis.buffer[j] = (delta + analysis.buffer[j]) * analysis.twiddles[i];
      }

      const size_t offsets[] = { kernelsize, kernelsize + (config.dftsize - 1) };

      for (size_t i = 1; i <= kernelsize; ++i)
      {
        analysis.buffer[offsets[0] - i] = std::conj(analysis.buffer[offsets[0] + i]);
        analysis.buffer[offsets[1] + i] = std::conj(analysis.buffer[offsets[1] - i]);
      }

      for (size_t i = 0; i < config.dftsize; ++i)
      {
        dft[i] = convolve(analysis.buffer.data() + i, analysis.window, analysis.weight);
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
        sdft(samples[i], &dfts[i * config.dftsize]);
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
        for (size_t i = 0; i < config.dftsize; ++i)
        {
          sample += dft[i].real() * (i % 2 ? -1 : +1);
        }
      }
      else
      {
        for (size_t i = 0; i < config.dftsize; ++i)
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
        samples[i] = isdft(&dfts[i * config.dftsize]);
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

    struct
    {
      double samplerate;
      size_t dftsize;
      std::vector<double> frequencies;
    }
    config;

    struct
    {
      Window window;
      F weight;
      struct { F first, last, all; } feedback;
      std::vector<std::complex<F>> buffer;
      std::vector<std::complex<F>> twiddles;
    }
    analysis;

    struct
    {
      double latency;
      F weight;
      std::vector<std::complex<F>> twiddles;
    }
    synthesis;

    static std::complex<F> convolve(const std::complex<F>* values, const Window window, const F weight)
    {
      switch (window)
      {
        case Window::Hann:
        {
          const std::complex<F> a = values[kernelsize] + values[kernelsize];
          const std::complex<F> b = values[kernelsize - 1] + values[kernelsize + 1];

          return (a - b) * weight * F(0.25);
        }
        case Window::Hamming:
        {
          const std::complex<F> a = F(0.54) * values[kernelsize];
          const std::complex<F> b = F(0.23) * (values[kernelsize - 1] + values[kernelsize + 1]);

          return (a - b) * weight;
        }
        case Window::Blackman:
        {
          const std::complex<F> a = F(0.42) * values[kernelsize];
          const std::complex<F> b = F(0.25) * (values[kernelsize - 1] + values[kernelsize + 1]);
          const std::complex<F> c = F(0.04) * (values[kernelsize - 2] + values[kernelsize + 2]);

          return (a - b + c) * weight;
        }
        default:
        {
          return values[kernelsize] * weight;
        }
      }
    }

  };
}
