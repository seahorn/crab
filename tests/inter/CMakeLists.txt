add_executable(inter test.cc)
target_link_libraries (inter ${TEST_LIBS})

install(TARGETS inter
  RUNTIME DESTINATION tests
  )
