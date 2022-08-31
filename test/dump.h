#pragma once

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

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

  void dump(const char* path, const sdft_double_complex_t* data, const size_t size)
  {
    FILE* file = fopen(path, "wb");
    fwrite(data, sizeof(sdft_double_complex_t), size, file);
    fclose(file);
  }

#endif

#if defined(_MSC_VER)
#pragma warning(pop)
#endif
