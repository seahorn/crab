
add_executable(test-1 test-1.cc)
target_link_libraries (test-1 ${CRAB_LIBS})

install(TARGETS test-1
  RUNTIME DESTINATION tests/domains/rationals
  )

add_executable(test-2 test-2.cc)
target_link_libraries (test-2 ${CRAB_LIBS})

install(TARGETS test-2
  RUNTIME DESTINATION tests/domains/rationals
  )


