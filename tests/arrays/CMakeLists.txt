add_executable(array_unit_tests unit-tests.cc)
target_link_libraries (array_unit_tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(array_tests tests.cc)
target_link_libraries (array_tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

