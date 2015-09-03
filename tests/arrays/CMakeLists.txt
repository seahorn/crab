add_executable(array_graph array_graph_tests.cc)
target_link_libraries (array_graph ${CRAB_LIBS})

add_executable(array_smashing array_smashing_tests.cc)
target_link_libraries (array_smashing ${CRAB_LIBS})
