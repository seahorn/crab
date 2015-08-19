add_executable(array_graph_tests array_graph_tests.cc)
target_link_libraries (array_graph_tests ${IKOS_CORE_LIBS})

add_executable(array_smashing_tests array_smashing_tests.cc)
target_link_libraries (array_smashing_tests ${IKOS_CORE_LIBS})
