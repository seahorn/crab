# GMP_USE_STATIC_LIBS   - Set to ON to force the use of static libraries.

## To switch between static and dynamic without resetting the cache
# if ("${GMP_USE_STATIC_LIBS}" STREQUAL "${GMP_USE_STATIC_LIBS_LAST}")
#   set(GMP_USE_STATIC_LIBS_CHANGED OFF)
# else ()
#   set(GMP_USE_STATIC_LIBS_CHANGED ON)
# endif()
# set(GMP_USE_STATIC_LIBS_LAST "${GMP_USE_STATIC_LIBS}")
# if (GMP_USE_STATIC_LIBS_CHANGED)
#   set(GMP_LIB "GMP_LIB-NOTFOUND")
# endif()

if (NOT GMP_FOUND)
  # Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
  if(GMP_USE_STATIC_LIBS )
    # save CMAKE_FIND_LIBRARY_SUFFIXES
    set(_GMP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()

  set(GMP_SEARCH_PATH "" CACHE PATH "Search path for gmp.")
  find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${GMP_SEARCH_PATH}/include)
  find_library(GMP_LIB NAMES gmp PATHS ${GMP_SEARCH_PATH}/lib)

  mark_as_advanced(GMP_SEARCH_PATH GMP_INCLUDE_DIR GMP_LIB)

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_INCLUDE_DIR GMP_LIB)

  if(GMP_USE_STATIC_LIBS )
    # restore CMAKE_FIND_LIBRARY_SUFFIXES
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_GMP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  
endif()
