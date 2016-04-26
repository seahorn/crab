add_executable(conc test1.cc)
target_link_libraries (conc ${CRAB_LIBS})

install(TARGETS conc
  RUNTIME DESTINATION tests/concurrency
  )


