add_executable(cg-tests cg-tests.cc)
target_link_libraries (cg-tests ${CRAB_LIBS})

install(TARGETS cg-tests
  RUNTIME DESTINATION tests/cg
  )


