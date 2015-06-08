add_executable(unit-tests unit-tests.cc)
target_link_libraries (unit-tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(tests tests.cc)
target_link_libraries (tests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

