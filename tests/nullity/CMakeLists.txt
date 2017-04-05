add_executable(null-tests null-tests.cc)
target_link_libraries (null-tests ${CRAB_LIBS})

install(TARGETS null-tests
  RUNTIME DESTINATION tests/domains
  )

