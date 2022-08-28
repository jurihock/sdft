#include <stdio.h>

#define SDFT_TD_FLOAT
#include <sdft/sdft.h>

#include <dump.h>
#include <wav.h>

int main(int argc, char* argv[])
{
  if (argc < 5)
  {
    return 1;
  }

  const size_t dftsize = atoi(argv[1]);
  const size_t hopsize = atoi(argv[2]);

  const char* ifile = argv[3];
  const char* ofile = argv[4];

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

  int progress = 0;

  for (size_t i = 0, j = 0; i < size; i+=hopsize, j++)
  {
    const double percent = (i / hopsize + 1.0) / (size / hopsize);

    if ((int)(percent * 10) != progress)
    {
        progress = (int)(percent * 10);
        printf("%i%%\n", progress * 10);
    }

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
