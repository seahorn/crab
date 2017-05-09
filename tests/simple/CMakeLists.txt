add_executable(test1 test1.cc)
target_link_libraries (test1 ${TEST_LIBS})

add_executable(test2 test2.cc)
target_link_libraries (test2 ${TEST_LIBS})

add_executable(test3 test3.cc)
target_link_libraries (test3 ${TEST_LIBS})

add_executable(test4 test4.cc)
target_link_libraries (test4 ${TEST_LIBS})

add_executable(test5 test5.cc)
target_link_libraries (test5 ${TEST_LIBS})

add_executable(unittests unittests.cc)
target_link_libraries (unittests ${TEST_LIBS})

# add_executable(sparsegraph sparsegraph.cc)
# target_link_libraries (sparsegraph ${TEST_LIBS})

install(TARGETS test1
  RUNTIME DESTINATION tests
  )

install(TARGETS test2
  RUNTIME DESTINATION tests
  )

install(TARGETS test3
  RUNTIME DESTINATION tests
  )

install(TARGETS test4
  RUNTIME DESTINATION tests
  )

install(TARGETS test5
  RUNTIME DESTINATION tests
  )


