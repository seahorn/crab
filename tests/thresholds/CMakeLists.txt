add_executable(thresholds test.cc)
target_link_libraries (thresholds ${CRAB_LIBS})

install(TARGETS thresholds
  RUNTIME DESTINATION tests/fixpo
  )


