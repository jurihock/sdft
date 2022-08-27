cmake_minimum_required(VERSION 3.1...3.18)

project(sdft)

if (MSVC)
  # add_compile_options(/W3 /WX)
else()
  # add_compile_options(-Wall -Werror)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/src/c/sdft/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/src/cpp/sdft/CMakeLists.txt")

include("${CMAKE_CURRENT_LIST_DIR}/lib/cpm.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/lib/dr_libs.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test/CMakeLists.txt")
