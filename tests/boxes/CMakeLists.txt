add_executable(ldd_test ldd_test.cc)
target_link_libraries (ldd_test ${CRAB_LIBS})

add_executable(boxes_test boxes_test.cc)
target_link_libraries (boxes_test ${CRAB_LIBS})
