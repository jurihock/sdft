#include <stdio.h>
#include <string.h>

#include <sdft/sdft.h>

#include <dump.h>
#include <wav.h>

sdft_window_t getwindow(const char* window)
{
  if (!strcmp(window, "hann"))
  {
    return sdft_window_hann;
  }

  if (!strcmp(window, "hamming"))
  {
    return sdft_window_hamming;
  }

  if (!strcmp(window, "blackman"))
  {
    return sdft_window_blackman;
  }

  return sdft_window_boxcar;
}

double getlatency(const char* latency)
{
  if (!strcmp(latency, "nan"))
  {
    return NAN;
  }

  return atoi(latency);
}

int main(int argc, char* argv[])
{
  if (argc < 8)
  {
    return 1;
  }

  const size_t dftsize = atoi(argv[1]);
  const size_t hopsize = atoi(argv[2]);
  const char* window = argv[3];
  const char* latency = argv[4];
  const char* srcfile = argv[5];
  const char* wavfile = argv[6];
  const char* dftfile = argv[7];

  sdft_t* sdft = sdft_alloc_custom(dftsize, getwindow(window), getlatency(latency));

  float* input;
  size_t size;
  size_t sr;

  if (!readwav(srcfile, &input, &size, &sr))
  {
    return 1;
  }

  printf("C\t%s %zu %zuHz\n", srcfile, size, sr);
  size = (size / hopsize) * hopsize;

  float* output = malloc(size * sizeof(float));
  sdft_double_complex_t* buffer = malloc(hopsize * dftsize * sizeof(sdft_double_complex_t));
  sdft_double_complex_t* dfts = malloc((size / hopsize) * dftsize * sizeof(sdft_double_complex_t));

  int progress = 0;

  for (size_t i = 0, j = 0; i < size; i+=hopsize, j++)
  {
    const double percent = (i / hopsize + 1.0) / (size / hopsize);

    if ((int)(percent * 10) != progress)
    {
        progress = (int)(percent * 10);
        printf("%i%%\n", progress * 10);
    }

    sdft_sdft_n(sdft, hopsize, input + i, buffer);
    sdft_isdft_n(sdft, hopsize, buffer, output + i);

    memcpy(dfts + j * dftsize, buffer, dftsize * sizeof(sdft_double_complex_t));
  }

  writewav(wavfile, output, size, sr);
  dump(dftfile, dfts, (size / hopsize) * dftsize);

  free(dfts);
  free(buffer);
  free(output);
  free(input);

  sdft_free(sdft);

  return 0;
}
