file(GLOB_RECURSE
  HDR "${CMAKE_CURRENT_LIST_DIR}/*.h")

# C

project(sdft-test-c)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE ${HDR} "${CMAKE_CURRENT_LIST_DIR}/test.c")

target_include_directories(${PROJECT_NAME}
  PRIVATE "${CMAKE_CURRENT_LIST_DIR}")

target_link_libraries(${PROJECT_NAME}
  PRIVATE sdft::c dr_libs)

# C++

project(sdft-test-cpp)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE ${HDR} "${CMAKE_CURRENT_LIST_DIR}/test.cpp")

target_include_directories(${PROJECT_NAME}
  PRIVATE "${CMAKE_CURRENT_LIST_DIR}")

target_link_libraries(${PROJECT_NAME}
  PRIVATE sdft::cpp dr_libs)
