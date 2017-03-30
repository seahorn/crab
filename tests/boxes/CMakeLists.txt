add_executable(ldd_test ldd_test.cc)
target_link_libraries (ldd_test ${CRAB_LIBS})

add_executable(boxes boxes_test.cc)
target_link_libraries (boxes ${CRAB_LIBS})

install(TARGETS boxes 
  RUNTIME DESTINATION tests/domains
  )
install(TARGETS ldd_test
  RUNTIME DESTINATION tests/domains
  )
