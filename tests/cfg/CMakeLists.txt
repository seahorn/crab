add_executable(cfg test.cc)
target_link_libraries (cfg ${TEST_LIBS})

install(TARGETS cfg
  RUNTIME DESTINATION tests
  )

