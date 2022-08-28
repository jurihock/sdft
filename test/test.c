#include <stdio.h>

#define SDFT_TD_FLOAT
#include <sdft/sdft.h>

#include <wav.h>

void dump(const char* path, const double complex* data, const size_t size)
{
  FILE* file = fopen(path, "wb");
  fwrite(data, sizeof(double complex), size, file);
  fclose(file);
}

int main()
{
  const size_t dftsize = 512;
  const size_t hopsize = 1000;

  const char* ifile = "test.wav";
  const char* ofile = "test.c.dfts";

  sdft_t* sdft = sdft_alloc(dftsize);

  float* samples;
  size_t size;
  double sr;

  if (!readwav(ifile, &samples, &size, &sr))
  {
    return 1;
  }

  printf("C\t%s %zu %gHz\n", ifile, size, sr);

  size = (size / hopsize) * hopsize;

  double complex* buffer = malloc(hopsize * dftsize * sizeof(double complex));
  double complex* dfts = malloc((size / hopsize) * dftsize * sizeof(double complex));

  for (size_t i = 0, j = 0; i < size; i+=hopsize, j++)
  {
    printf("%zu/%zu\n", i / hopsize + 1, size / hopsize);

    sdft_sdft_n(sdft, hopsize, samples + i, buffer);

    memcpy(dfts + j * dftsize, buffer, dftsize * sizeof(double complex));
  }

  dump(ofile, dfts, (size / hopsize) * dftsize);

  free(dfts);
  free(buffer);
  free(samples);
  sdft_free(sdft);

  return 0;
}
