add_executable(checkers test.cc)
target_link_libraries (checkers ${TEST_LIBS})

install(TARGETS checkers
  RUNTIME DESTINATION tests
  )
