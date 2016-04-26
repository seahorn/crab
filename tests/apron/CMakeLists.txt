add_executable(apron-domains apron_test.cc)
target_link_libraries (apron-domains ${CRAB_LIBS})

install(TARGETS apron-domains
  RUNTIME DESTINATION tests/domains
  )


