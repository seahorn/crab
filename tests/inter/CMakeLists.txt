add_executable(inter-tests inter-tests.cc)
target_link_libraries (inter-tests ${CRAB_LIBS})

install(TARGETS inter-tests
  RUNTIME DESTINATION tests/inter
  )
