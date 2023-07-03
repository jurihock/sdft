#include <sdft/sdft.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <iostream>
#include <vector>

using sdft::SDFT;

int main()
{
  std::cout << "CPP;\tSDFT;\tISDFT" << std::endl;

  const auto samplerate = 44100;
  const auto dftsize = 1000;

  const auto ta0 = std::chrono::high_resolution_clock::now();
  SDFT<double, double> sdft(dftsize);
  const auto tb0 = std::chrono::high_resolution_clock::now();
  const auto e0 = std::chrono::duration_cast<std::chrono::microseconds>(tb0 - ta0).count();

  std::cout << "0" << ";\t" << e0 << ";\t" << e0 << std::endl;

  const auto n = 1 * samplerate;
  const auto m = sdft.size();

  std::vector<double> x(n);
  std::vector<std::complex<double>> y(n * m);

  const auto runs = 10;

  for (auto run = 1; run <= runs; ++run)
  {
    std::fill(x.begin(), x.end(), 0.0);
    std::fill(y.begin(), y.end(), 0.0);

    const auto ta1 = std::chrono::high_resolution_clock::now();
    sdft.sdft(n, x.data(), y.data());
    const auto tb1 = std::chrono::high_resolution_clock::now();
    const auto e1 = std::chrono::duration_cast<std::chrono::microseconds>(tb1 - ta1).count();

    const auto ta2 = std::chrono::high_resolution_clock::now();
    sdft.isdft(n, y.data(), x.data());
    const auto tb2 = std::chrono::high_resolution_clock::now();
    const auto e2 = std::chrono::duration_cast<std::chrono::microseconds>(tb2 - ta2).count();

    std::cout << run << ";\t" << e1 << ";\t" << e2 << std::endl;
  }

  return 0;
}
