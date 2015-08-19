add_executable(simple_test1 prog-1.cc)
target_link_libraries (simple_test1 ${IKOS_CORE_LIBS})

add_executable(simple_test2 prog-2.cc)
target_link_libraries (simple_test2 ${IKOS_CORE_LIBS})

add_executable(simple_test3 prog-3.cc)
target_link_libraries (simple_test3 ${IKOS_CORE_LIBS})

add_executable(simple_test4 prog-4.cc)
target_link_libraries (simple_test4 ${IKOS_CORE_LIBS})

add_executable(unittests unittests.cc)
target_link_libraries (unittests ${IKOS_CORE_LIBS})

