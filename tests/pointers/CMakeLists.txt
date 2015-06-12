add_executable(pointer_tests tests.cc)
target_link_libraries (pointer_tests ${GMPXX_LIB} ${GMP_LIB} pthread)

