# Find Cudd

if (NOT CUDD_ROOT)
  set(CUDD_ROOT "" CACHE PATH "Root of Cudd compiled source tree.")
endif()

find_path(CUDD_INCLUDE_DIR NAMES cudd.h PATHS ${CUDD_ROOT}/include/ldd)
#find_path(CUDD_INT_INCLUDE_DIR NAMES cuddInt.h PATHS ${CUDD_ROOT}/cudd)
find_path(CUDD_INT_INCLUDE_DIR NAMES cuddInt.h PATHS ${CUDD_ROOT}/include/ldd)

set (CUDD_INCLUDE_DIR ${CUDD_INCLUDE_DIR} ${CUDD_INT_INCLUDE_DIR})

mark_as_advanced (CUDD_INCLUDE_DIR CUDD_INT_INCLUDE_DIR)

find_library(CUDD_CUDD_LIBRARY   NAMES cudd   PATHS  ${CUDD_ROOT}/lib)
find_library(CUDD_EPD_LIBRARY   NAMES epd   PATHS  ${CUDD_ROOT}/lib)
find_library(CUDD_ST_LIBRARY   NAMES st   PATHS  ${CUDD_ROOT}/lib)
## Hack: In some linux systems, the library /lib64/libutil.so exists
## but it has nothing to do with cudd. Even if we use NO_DEFAULT_PATH
## when CUDD_ROOT is not defined CUDD_UTIL_LIBRARY can be set to the
## wrong library.
if (CUDD_CUDD_LIBRARY)
   find_library(CUDD_UTIL_LIBRARY util PATHS ${CUDD_ROOT}/lib NO_DEFAULT_PATH)
endif()
find_library(CUDD_MTR_LIBRARY NAMES mtr PATHS ${CUDD_ROOT}/lib)

mark_as_advanced (CUDD_CUDD_LIBRARY CUDD_DDDMP_LIBRARY CUDD_EPD_LIBRARY
  CUDD_ST_LIBRARY CUDD_UTIL_LIBRARY CUDD_MTR_LIBRARY)

set(CUDD_LIBRARY 
  ${CUDD_CUDD_LIBRARY} ${CUDD_ST_LIBRARY} ${CUDD_UTIL_LIBRARY}
  ${CUDD_MTR_LIBRARY}  ${CUDD_EPD_LIBRARY} ${CUDD_DDDMP_LIBRARY} )
mark_as_advanced (CUDD_LIBRARY)

include (CheckTypeSize)
check_type_size (long CMAKE_SIZEOF_LONG)
#message (STATUS "sizeof (long): ${CMAKE_SIZEOF_LONG}")

set (CUDD_CXXFLAGS "-DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=${CMAKE_SIZEOF_VOID_P} -DSIZEOF_LONG=${CMAKE_SIZEOF_LONG}")
mark_as_advanced (CMAKE_SIZEOF_LONG CUDD_CXXFLAGS)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cudd
  REQUIRED_VARS CUDD_LIBRARY CUDD_INCLUDE_DIR CUDD_CXXFLAGS)

