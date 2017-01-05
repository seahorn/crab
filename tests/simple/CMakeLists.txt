add_executable(test1 test1.cc)
target_link_libraries (test1 ${CRAB_LIBS})

add_executable(test2 test2.cc)
target_link_libraries (test2 ${CRAB_LIBS})

add_executable(test3 test3.cc)
target_link_libraries (test3 ${CRAB_LIBS})

add_executable(test4 test4.cc)
target_link_libraries (test4 ${CRAB_LIBS})

add_executable(test5 test5.cc)
target_link_libraries (test5 ${CRAB_LIBS})

add_executable(test6 test6.cc)
target_link_libraries (test6 ${CRAB_LIBS})

add_executable(unittests unittests.cc)
target_link_libraries (unittests ${CRAB_LIBS})

# add_executable(sparsegraph sparsegraph.cc)
# target_link_libraries (sparsegraph ${CRAB_LIBS})

install(TARGETS test1
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS test2
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS test3
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS test4
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS test5
  RUNTIME DESTINATION tests/domains
  )

install(TARGETS test6
  RUNTIME DESTINATION tests/domains
  )


