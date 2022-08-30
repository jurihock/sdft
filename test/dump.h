#pragma once

#if defined(__cplusplus)

  #include <cstdio>

  void dump(const char* path, const std::complex<double>* data, const size_t size)
  {
    FILE* file = fopen(path, "wb");
    fwrite(data, sizeof(std::complex<double>), size, file);
    fclose(file);
  }

#else

  #include <stdio.h>

  void dump(const char* path, const sdft_fdx_t* data, const size_t size)
  {
    FILE* file = fopen(path, "wb");
    fwrite(data, sizeof(sdft_fdx_t), size, file);
    fclose(file);
  }

#endif
