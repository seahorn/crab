add_executable(disintervals-tests disintervals-tests.cc)
target_link_libraries (disintervals-tests ${CRAB_LIBS})

install(TARGETS disintervals-tests
  RUNTIME DESTINATION tests/domains
  )


