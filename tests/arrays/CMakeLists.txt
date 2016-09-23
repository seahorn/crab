add_executable(array_smashing array_smashing_tests.cc)
target_link_libraries (array_smashing ${CRAB_LIBS})

add_executable(array_sparse_graph array_sparse_graph_tests.cc)
target_link_libraries (array_sparse_graph ${CRAB_LIBS})

install(TARGETS array_smashing
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS array_sparse_graph
  RUNTIME DESTINATION tests/domains
  )

