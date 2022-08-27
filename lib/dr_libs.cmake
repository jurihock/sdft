# https://github.com/mackron/dr_libs

CPMAddPackage(
  NAME dr_libs
  VERSION 0.13.6
  GIT_TAG 15f37e3ab01654c1a3bc98cff2a9ca64e8296fa9
  GITHUB_REPOSITORY mackron/dr_libs
  DOWNLOAD_ONLY YES)

if(dr_libs_ADDED)

  add_library(dr_libs INTERFACE)

  target_include_directories(dr_libs
    INTERFACE "${dr_libs_SOURCE_DIR}")

  target_compile_definitions(dr_libs
    INTERFACE -DDR_WAV_IMPLEMENTATION)

endif()
