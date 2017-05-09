# add_executable(ldd_test ldd_test.cc)
# target_link_libraries (ldd_test ${TEST_LIBS})

add_executable(boxes boxes_test.cc)
target_link_libraries (boxes ${TEST_LIBS})

install(TARGETS boxes
  RUNTIME DESTINATION tests
  )
