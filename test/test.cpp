#include <cstdio>

#include <sdft/sdft.h>

#include <dump.h>
#include <wav.h>

int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    return 1;
  }

  const size_t dftsize = atoi(argv[1]);
  const size_t hopsize = atoi(argv[2]);

  const char* ifile = argv[3];
  const char* ofile = argv[4];

  SDFT<float> sdft(dftsize);

  float* samples;
  size_t size;
  size_t sr;

  if (!readwav(ifile, &samples, &size, &sr))
  {
    return 1;
  }

  printf("C++\t%s %zu %zuHz\n", ifile, size, sr);
  size = (size / hopsize) * hopsize;

  std::complex<double>* buffer = new std::complex<double>[hopsize * dftsize];
  std::complex<double>* dfts = new std::complex<double>[(size / hopsize) * dftsize];

  int progress = 0;

  for (size_t i = 0, j = 0; i < size; i+=hopsize, j++)
  {
    const double percent = (i / hopsize + 1.0) / (size / hopsize);

    if ((int)(percent * 10) != progress)
    {
        progress = (int)(percent * 10);
        printf("%i%%\n", progress * 10);
    }

    sdft.sdft(hopsize, samples + i, buffer);
    memcpy(dfts + j * dftsize, buffer, dftsize * sizeof(std::complex<double>));
  }

  dump(ofile, dfts, (size / hopsize) * dftsize);

  delete[] dfts;
  delete[] buffer;
  delete[] samples;

  return 0;
}
