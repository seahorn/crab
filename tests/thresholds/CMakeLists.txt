add_executable(thresholds test.cc)
target_link_libraries (thresholds ${TEST_LIBS})

install(TARGETS thresholds
  RUNTIME DESTINATION tests
  )


