add_executable(simple_test1 prog-1.cc)
target_link_libraries (simple_test1 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(simple_test2 prog-2.cc)
target_link_libraries (simple_test2 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(simple_test3 prog-3.cc)
target_link_libraries (simple_test3 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(simple_test4 prog-4.cc)
target_link_libraries (simple_test4 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(unittests unittests.cc)
target_link_libraries (unittests ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

