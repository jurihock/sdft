#include <stdio.h>

#define SDFT_TD_FLOAT
#include <sdft/sdft.h>

#include <wav.h>

int main()
{
  const size_t dftsize = 512;
  const size_t hopsize = 1000;
  const char* file = "test.wav";

  sdft_t* sdft = sdft_alloc(dftsize);

  float* data;
  size_t size;
  double sr;

  if (!readwav(file, &data, &size, &sr))
  {
    return 1;
  }

  printf("C\t%s %zu %gHz\n", file, size, sr);

  size = (size / hopsize) * hopsize;

  const size_t nsamples = size / hopsize;

  double complex* dfts = malloc(nsamples * dftsize * sizeof(double complex));

  for (size_t i = 0; i < size; i+=hopsize)
  {
    printf("%zu/%zu\n", i / hopsize + 1, size / hopsize);

    float* samples = data + i;

    sdft_sdft_n(sdft, nsamples, samples, dfts);
  }

  free(dfts);
  free(data);
  sdft_free(sdft);

  return 0;
}
