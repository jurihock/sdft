cmake_minimum_required(VERSION 3.12...3.18)

project(sdft-test)

if(MSVC)
  # add_compile_options(/W3 /WX)
else()
  # add_compile_options(-Wall -Werror)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../cpp/lib/cpm.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../cpp/lib/dr.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/../c/src/sdft/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/../cpp/src/sdft/CMakeLists.txt")

foreach(XXX c cpp)

  project(sdft-test-${XXX})

  add_executable(${PROJECT_NAME})

  file(GLOB
    HDR "${CMAKE_CURRENT_LIST_DIR}/*.h")

  file(GLOB
    SRC "${CMAKE_CURRENT_LIST_DIR}/*.${XXX}")

  target_sources(${PROJECT_NAME}
    PRIVATE ${HDR} "${SRC}")

  target_include_directories(${PROJECT_NAME}
    PRIVATE "${CMAKE_CURRENT_LIST_DIR}")

  target_link_libraries(${PROJECT_NAME}
    PRIVATE sdft::${XXX} dr)

  if (MSVC)
    target_compile_options(${PROJECT_NAME}
      PRIVATE /fp:fast)
  else()
    target_compile_options(${PROJECT_NAME}
      PRIVATE -ffast-math)
  endif()

endforeach()
