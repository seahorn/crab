add_library(Debug ${CRAB_LIBS_TYPE} 
  debug.cpp)

install(TARGETS Debug
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
