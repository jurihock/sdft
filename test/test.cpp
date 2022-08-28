#include <cstdio>

#include <sdft/sdft.h>

#include <dump.h>
#include <wav.h>

int main()
{
  const size_t dftsize = 512;
  const size_t hopsize = 1000;

  const char* ifile = "test.wav";
  const char* ofile = "test.cpp.dfts";

  SDFT<float> sdft(dftsize);

  float* samples;
  size_t size;
  double sr;

  if (!readwav(ifile, &samples, &size, &sr))
  {
    return 1;
  }

  printf("C++\t%s %zu %gHz\n", ifile, size, sr);
  size = (size / hopsize) * hopsize;

  std::complex<double>* buffer = new std::complex<double>[hopsize * dftsize];
  std::complex<double>* dfts = new std::complex<double>[(size / hopsize) * dftsize];

  for (size_t i = 0, j = 0; i < size; i+=hopsize, j++)
  {
    printf("%zu/%zu\n", i / hopsize + 1, size / hopsize);
    sdft.sdft(hopsize, samples + i, buffer);
    memcpy(dfts + j * dftsize, buffer, dftsize * sizeof(std::complex<double>));
  }

  dump(ofile, dfts, (size / hopsize) * dftsize);

  delete[] dfts;
  delete[] buffer;
  delete[] samples;

  return 0;
}
