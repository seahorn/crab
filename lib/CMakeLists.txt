add_library(Crab ${CRAB_LIBS_TYPE} 
  debug.cpp
  stats.cpp
  os.cpp
  term.cpp
  )

install(TARGETS Crab
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

