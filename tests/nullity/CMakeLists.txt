add_executable(nullity test.cc)
target_link_libraries (nullity ${CRAB_LIBS})

install(TARGETS nullity
  RUNTIME DESTINATION tests/domains
  )

