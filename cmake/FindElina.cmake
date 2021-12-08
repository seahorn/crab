#Find Elina library
if (NOT ELINA_FOUND)
  
  # save CMAKE_FIND_LIBRARY_SUFFIXES
  set(_ELINA_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

  ### Elina does not produce static libraries  
  #if(NOT BUILD_CRAB_LIBS_SHARED)
  #  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  #else()
  ### Elina dynamic libraries have .so extension even if they are built
  ### on mac osx.
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib)      
  #endif()
  
  set (ELINA_ROOT "" CACHE PATH "Root of Elina install.")
  
  find_package (Gmp QUIET)
  find_package (MPFR QUIET)
  
  find_path(ELINA_INCLUDE_DIR NAMES elina_abstract0.h PATHS ${ELINA_ROOT}/include)
  
  find_library(Elina_Poly_Lib NAMES optpoly PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Elina_Opt_Oct_Lib NAMES optoct PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Elina_Linearize_Lib NAMES elinalinearize PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Elina_Aux_Lib NAMES elinaux PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Elina_Opt_Zones_Lib NAMES optzones PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Elina_Partitions_Lib NAMES partitions PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
  
  set(ELINA_LIBRARY ${Elina_Poly_Lib} 
    ${Elina_Opt_Oct_Lib} ${Elina_Linearize_Lib} ${Elina_Aux_Lib}
    ${Elina_Opt_Zones_Lib} ${Elina_Partitions_Lib})
  
  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (Elina
    REQUIRED_VARS ELINA_INCLUDE_DIR ELINA_LIBRARY GMP_FOUND MPFR_FOUND)
  
  set (ELINA_INCLUDE_DIR ${ELINA_INCLUDE_DIR} ${MPFR_INC_DIR})
  set (ELINA_LIBRARY ${ELINA_LIBRARY} ${MPFR_LIB})
  
  mark_as_advanced(ELINA_LIBRARY ELINA_INCLUDE_DIR)
  
  # restore CMAKE_FIND_LIBRARY_SUFFIXES
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ELINA_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
   
endif ()
