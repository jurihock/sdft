cmake_minimum_required(VERSION 3.12...3.18)

project(sdft)

if(MSVC)
  # add_compile_options(/W3 /WX)
else()
  # add_compile_options(-Wall -Werror)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/src/sdft/CMakeLists.txt")
