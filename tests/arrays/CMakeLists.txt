add_executable(array_graph array_graph_tests.cc)
target_link_libraries (array_graph ${TEST_LIBS})

add_executable(array_smashing array_smashing_tests.cc)
target_link_libraries (array_smashing ${TEST_LIBS})

install(TARGETS array_graph
  RUNTIME DESTINATION tests
  )

install(TARGETS array_smashing
  RUNTIME DESTINATION tests
  )
