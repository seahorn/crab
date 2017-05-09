add_executable(callgraph test.cc)
target_link_libraries (callgraph ${TEST_LIBS})

install(TARGETS callgraph
  RUNTIME DESTINATION tests
  )


