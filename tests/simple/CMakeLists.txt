add_executable(test1 test1.cc)
target_link_libraries (test1 ${CRAB_LIBS})

add_executable(test2 test2.cc)
target_link_libraries (test2 ${CRAB_LIBS})

add_executable(test3 test3.cc)
target_link_libraries (test3 ${CRAB_LIBS})

add_executable(test4 test4.cc)
target_link_libraries (test4 ${CRAB_LIBS})

add_executable(unittests unittests.cc)
target_link_libraries (unittests ${CRAB_LIBS})

