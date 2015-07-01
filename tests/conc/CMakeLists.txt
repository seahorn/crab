add_executable(conc_test1 test1.cc)
target_link_libraries (conc_test1 ${GMPXX_LIB} ${GMP_LIB} dbm pthread)


