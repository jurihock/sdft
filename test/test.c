#include <stdio.h>

#define SDFT_TD_FLOAT
#include <sdft/sdft.h>

int main()
{
  SDFT* sdft = sdft_alloc(1024);

  sdft_free(sdft);

  return 0;
}
