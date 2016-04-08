add_library(Stats ${CRAB_LIBS_TYPE} 
  stats.cpp)

install(TARGETS Stats
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
