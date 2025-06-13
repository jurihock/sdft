# https://github.com/mackron/dr_libs

CPMAddPackage(
  NAME dr
  VERSION 2025-05-30
  GIT_TAG 92844c5c07f05b21855e37b482ed0b3143256cf6
  GITHUB_REPOSITORY mackron/dr_libs
  DOWNLOAD_ONLY YES)

if(dr_ADDED)

  add_library(dr INTERFACE)

  target_include_directories(dr
    INTERFACE "${dr_SOURCE_DIR}")

  target_compile_definitions(dr
    INTERFACE -DDR_WAV_IMPLEMENTATION)

endif()
