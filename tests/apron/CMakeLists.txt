add_executable(apron-domains apron_test.cc)
target_link_libraries (apron-domains ${TEST_LIBS})

install(TARGETS apron-domains
  RUNTIME DESTINATION tests
  )


