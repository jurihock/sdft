#include <iostream>

#include <sdft/sdft.h>

#include <wav.h>

int main()
{
  const size_t dftsize = 512;
  const size_t hopsize = 1000;
  const char* file = "test.wav";

  SDFT<float> sdft(dftsize);

  float* data;
  size_t size;
  double sr;

  if (!readwav(file, &data, &size, &sr))
  {
    return 1;
  }

  std::cout << "C++\t" << file << " " << size << " " << sr << "Hz" << std::endl;

  size = (size / hopsize) * hopsize;

  const size_t nsamples = size / hopsize;
  std::complex<double>* dfts = new std::complex<double>[nsamples * dftsize];

  for (size_t i = 0; i < size; i+=hopsize)
  {
    std::cout << i / hopsize + 1 << "/" << size / hopsize << std::endl;

    float* samples = data + i;

    sdft.sdft(nsamples, samples, dfts);
  }

  delete[] dfts;
  delete[] data;

  return 0;
}
