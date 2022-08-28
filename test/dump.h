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

  void dump(const char* path, const double complex* data, const size_t size)
  {
    FILE* file = fopen(path, "wb");
    fwrite(data, sizeof(double complex), size, file);
    fclose(file);
  }

#endif
