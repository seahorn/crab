add_library(term ${CRAB_LIBS_TYPE}  
  term.cpp)


install(TARGETS term
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
