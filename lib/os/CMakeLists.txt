add_library(OS ${CRAB_LIBS_TYPE} 
  os.cpp)

install(TARGETS OS
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
