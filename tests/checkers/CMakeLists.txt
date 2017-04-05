add_executable(checker-tests checker-tests.cc)
target_link_libraries (checker-tests ${CRAB_LIBS})

install(TARGETS checker-tests
  RUNTIME DESTINATION tests/checkers
  )
