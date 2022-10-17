# https://github.com/soblin/matplotlibcpp17

CPMAddPackage(
  NAME matplotlibcpp
  VERSION 1.0.4
  GIT_TAG bb9e692c1b48b83069e2bb0e17354098b841304f
  GITHUB_REPOSITORY soblin/matplotlibcpp17
  DOWNLOAD_ONLY YES)

if(matplotlibcpp_ADDED)

  add_library(matplotlibcpp INTERFACE)

  target_include_directories(matplotlibcpp
    INTERFACE "${matplotlibcpp_SOURCE_DIR}/include")

  target_compile_features(matplotlibcpp
    INTERFACE cxx_std_17)

endif()
