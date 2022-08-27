project(sdft-lib-cpp)

add_library(${PROJECT_NAME} INTERFACE)
add_library(sdft::cpp ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  INTERFACE "${CMAKE_CURRENT_LIST_DIR}/sdft.h")

target_include_directories(${PROJECT_NAME}
  INTERFACE "${CMAKE_CURRENT_LIST_DIR}/..")

target_compile_features(${PROJECT_NAME}
  INTERFACE cxx_std_11)
