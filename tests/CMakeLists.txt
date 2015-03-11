add_executable(prog-1 prog-1.cc)
target_link_libraries (prog-1 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(prog-2 prog-2.cc)
target_link_libraries (prog-2 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(prog-3 prog-3.cc)
target_link_libraries (prog-3 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)

add_executable(unif-1 unif-1.cc)
target_link_libraries (unif-1 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)
