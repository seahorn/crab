add_library(Term ${CRAB_LIBS_TYPE}  
  term.cpp)


install(TARGETS Term
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
