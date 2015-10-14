add_library(boxes STATIC boxes.cpp)


install(TARGETS boxes
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
