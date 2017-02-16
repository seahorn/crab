add_executable(apron-domains apron_test.cc)
target_link_libraries (apron-domains ${CRAB_LIBS})

add_executable(apron-bug apron_bug.cc)
target_link_libraries (apron-bug ${CRAB_LIBS})

install(TARGETS apron-domains apron-bug
  RUNTIME DESTINATION tests/domains
  )


