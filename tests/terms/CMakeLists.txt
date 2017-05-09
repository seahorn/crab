# add_executable(term-0 terms-0.cc)
# target_link_libraries (term-0 ${TEST_LIBS})

# add_executable(term-0a terms-0a.cc)
# target_link_libraries (term-0a ${TEST_LIBS})

add_executable(term-1 terms-1.cc)
target_link_libraries (term-1 ${TEST_LIBS})

add_executable(term-2 terms-2.cc)
target_link_libraries (term-2 ${TEST_LIBS})

add_executable(term-3 terms-3.cc)
target_link_libraries (term-3 ${TEST_LIBS})

add_executable(term-4 terms-4.cc)
target_link_libraries (term-4 ${TEST_LIBS})

install(TARGETS term-1
  RUNTIME DESTINATION tests
  )

install(TARGETS term-2
  RUNTIME DESTINATION tests
  )

install(TARGETS term-3
  RUNTIME DESTINATION tests
  )

install(TARGETS term-4
  RUNTIME DESTINATION tests
  )

