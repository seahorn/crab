#Find Apron library
if (NOT APRON_FOUND)
  
  # save CMAKE_FIND_LIBRARY_SUFFIXES
  set(_APRON_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

  #### FIXME: we cannot use apron dynamic libraries because crab
  #### library refers to the apron libraries located in the build
  #### directory instead of the ones where they will be installed.
  #if(NOT BUILD_CRAB_LIBS_SHARED)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  #else()
  #  set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib)
  #endif()
  
  set (APRON_ROOT "" CACHE PATH "Root of Apron install.")
  
  find_package (Gmp QUIET)
  find_package (MPFR QUIET)
  
  find_path(APRON_INCLUDE_DIR NAMES ap_abstract0.h PATHS ${APRON_ROOT}/include)
  
  find_library(Apron_Polka_Lib NAMES polkaMPQ PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  ## octD is faster than octMPQ
  find_library(Apron_Oct_Lib NAMES octD PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  # find_library(Apron_Oct_Lib NAMES octMPQ PATHS ${APRON_ROOT}/lib)
  find_library(Apron_Opt_Oct_Lib NAMES optoct PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Apron_Opt_Oct_utils_Lib NAMES linkedlistapi PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Apron_Apron_Lib NAMES apron PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  find_library(Apron_Box_Lib NAMES boxMPQ PATHS ${APRON_ROOT}/lib NO_DEFAULT_PATH)
  
  set(APRON_LIBRARY ${Apron_Box_Lib} 
    ${Apron_Polka_Lib} ${Apron_Oct_Lib} 
    ${Apron_Opt_Oct_Lib} ${Apron_Opt_Oct_utils_Lib}
    ${Apron_Apron_Lib})
  
  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (Apron
    REQUIRED_VARS APRON_INCLUDE_DIR APRON_LIBRARY GMP_FOUND MPFR_FOUND)
  
  set (APRON_INCLUDE_DIR ${APRON_INCLUDE_DIR} ${MPFR_INC_DIR})
  set (APRON_LIBRARY ${APRON_LIBRARY} ${MPFR_LIB})
  
  mark_as_advanced(APRON_LIBRARY APRON_INCLUDE_DIR 
    Apron_Apron_Lib Apron_Box_Lib Apron_Oct_Lib Apron_Polka_Lib)
  
  # restore CMAKE_FIND_LIBRARY_SUFFIXES
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_APRON_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
  
endif ()
