add_library(term STATIC term.cpp)


install(TARGETS term
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
