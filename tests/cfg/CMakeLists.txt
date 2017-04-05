add_executable(cfg-tests cfg-tests.cc)
target_link_libraries (cfg-tests ${CRAB_LIBS})

install(TARGETS cfg-tests
  RUNTIME DESTINATION tests/cfg
  )

