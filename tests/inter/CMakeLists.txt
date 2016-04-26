add_executable(inter test.cc)
target_link_libraries (inter ${CRAB_LIBS})

install(TARGETS inter
  RUNTIME DESTINATION tests/inter
  )
