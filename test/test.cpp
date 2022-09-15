#include <cstdio>
#include <cstring>

#include <sdft/sdft.h>

#include <dump.h>
#include <wav.h>

SDFT<>::Window getwindow(const char* window)
{
  if (!strcmp(window, "hann"))
  {
    return SDFT<>::Window::Hann;
  }

  if (!strcmp(window, "hamming"))
  {
    return SDFT<>::Window::Hamming;
  }

  return SDFT<>::Window::Boxcar;
}

int main(int argc, char *argv[])
{
  if (argc < 7)
  {
    return 1;
  }

  const size_t dftsize = atoi(argv[1]);
  const size_t hopsize = atoi(argv[2]);
  const char* window = argv[3];

  const char* srcfile = argv[4];
  const char* wavfile = argv[5];
  const char* dftfile = argv[6];

  SDFT<> sdft(dftsize, getwindow(window));

  float* input;
  size_t size;
  size_t sr;

  if (!readwav(srcfile, &input, &size, &sr))
  {
    return 1;
  }

  printf("C++\t%s %zu %zuHz\n", srcfile, size, sr);
  size = (size / hopsize) * hopsize;

  float* output = new float[size];
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

    sdft.sdft(hopsize, input + i, buffer);
    sdft.isdft(hopsize, buffer, output + i);

    memcpy(dfts + j * dftsize, buffer, dftsize * sizeof(std::complex<double>));
  }

  writewav(wavfile, output, size, sr);
  dump(dftfile, dfts, (size / hopsize) * dftsize);

  delete[] dfts;
  delete[] buffer;
  delete[] output;
  delete[] input;

  return 0;
}
