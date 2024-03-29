if(Python3_FOUND AND Python3_NumPy_FOUND)

  project(sdft-example-analysis)

  add_executable(${PROJECT_NAME})

  target_sources(${PROJECT_NAME}
    PRIVATE "${CMAKE_CURRENT_LIST_DIR}/analysis.cpp")

  target_link_libraries(${PROJECT_NAME}
    PRIVATE sdft::cpp matplotlibcpp numcpp pybind python)

  if (MSVC)
    target_compile_options(${PROJECT_NAME}
      PRIVATE /fp:fast)
  else()
    target_compile_options(${PROJECT_NAME}
      PRIVATE -ffast-math)
  endif()

endif()
