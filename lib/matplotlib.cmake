# https://github.com/lava/matplotlib-cpp

CPMAddPackage(
  NAME matplotlib
  VERSION 2021.04.23
  GIT_TAG ef0383f1315d32e0156335e10b82e90b334f6d9f
  GITHUB_REPOSITORY lava/matplotlib-cpp
  DOWNLOAD_ONLY YES)

if(matplotlib_ADDED)

  file(
    COPY "${CMAKE_CURRENT_LIST_DIR}/matplotlib.patch"
    DESTINATION "${matplotlib_SOURCE_DIR}")

  execute_process(
    COMMAND git apply --quiet matplotlib.patch
    WORKING_DIRECTORY "${matplotlib_SOURCE_DIR}")

  add_library(matplotlib INTERFACE)

  target_include_directories(matplotlib
    INTERFACE "${matplotlib_SOURCE_DIR}")

  target_compile_features(matplotlib
    INTERFACE cxx_std_17)

endif()
