add_executable(thresholds-tests thresholds-tests.cc)
target_link_libraries (thresholds-tests ${CRAB_LIBS})

install(TARGETS thresholds-tests
  RUNTIME DESTINATION tests/fixpo
  )


