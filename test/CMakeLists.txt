# C

project(sdft-test-c)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.c")

target_link_libraries(${PROJECT_NAME}
  PRIVATE sdft::c)

# C++

project(sdft-test-cpp)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.cpp")

target_link_libraries(${PROJECT_NAME}
  PRIVATE sdft::cpp)
