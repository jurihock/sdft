foreach(XXX c cpp)

  project(sdft-test-${XXX})

  add_executable(${PROJECT_NAME})

  file(GLOB_RECURSE
    HDR "${CMAKE_CURRENT_LIST_DIR}/*.h")

  file(GLOB_RECURSE
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
