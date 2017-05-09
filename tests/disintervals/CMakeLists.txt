add_executable(dis-int test.cc)
target_link_libraries (dis-int ${TEST_LIBS})

install(TARGETS dis-int
  RUNTIME DESTINATION tests
  )


