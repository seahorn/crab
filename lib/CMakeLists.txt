add_library(Crab ${CRAB_LIBS_TYPE} 
  debug.cpp
  stats.cpp
  os.cpp
  boolean.cpp
  ##numbers 
  bignums.cpp
  safeint.cpp
  wrapint.cpp
  ##variable factory
  varname_factory.cpp
  ##types
  tag.cpp
  )

target_link_libraries(Crab ${GMP_LIB})

install(TARGETS Crab
  LIBRARY DESTINATION crab/lib
  ARCHIVE DESTINATION crab/lib)

if (NOT TopLevel)
  ## This is ugly fix so that seahorn and crab-llvm can find Crab if
  ## dynamically linked. It has to do with how CMAKE_INSTALL_RPATH is
  ## handled by these tools.
  install(TARGETS Crab
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()  



