add_executable(array_graph_tests array_graph_tests.cc)
target_link_libraries (array_graph_tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(array_smashing_tests array_smashing_tests.cc)
target_link_libraries (array_smashing_tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)
