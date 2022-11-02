# Find PPLite headers and libraries.
# Author: Enea Zaffanella (enea.zaffanella@unipr.it)

if (NOT PPLITE_FOUND)
  set(PPLITE_ROOT "" CACHE PATH "Path to PPLite install directory")

  set(PPLITE_INCLUDE_SEARCH_DIRS "")
  set(PPLITE_LIB_SEARCH_DIRS "")

  if (PPLITE_ROOT)
    list(APPEND PPLITE_INCLUDE_SEARCH_DIRS "${PPLITE_ROOT}/include/")
    list(APPEND PPLITE_LIB_SEARCH_DIRS "${PPLITE_ROOT}/lib")
  endif()

  find_package(GMP QUIET)
  find_package(GMPXX QUIET)
  find_package(FLINT QUIET)
  cmake_policy(SET CMP0074 NEW)
  find_package(APRON QUIET)

  find_path(PPLITE_INCLUDE_DIR
    NAMES pplite/pplite.hh
    HINTS ${PPLITE_INCLUDE_SEARCH_DIRS}
    DOC "Path to PPLite include directory"
  )

  find_library(PPLITE_LIB
    NAMES pplite
    HINTS ${PPLITE_LIB_SEARCH_DIRS}
    DOC "Path to PPLite library"
  )

  find_path(AP_PPLITE_INCLUDE_DIR
    NAMES ap_pplite.h
    HINTS ${PPLITE_INCLUDE_SEARCH_DIRS}
    DOC "Path to Apron PPLite's wrapper include directory"
  )

  find_library(AP_PPLITE_LIB
    NAMES ap_pplite
    HINTS ${PPLITE_LIB_SEARCH_DIRS}
    DOC "Path to Apron PPLite's wrapper library ap_pplite"
  )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PPLITE
    REQUIRED_VARS
      PPLITE_INCLUDE_DIR
      PPLITE_LIB
      AP_PPLITE_INCLUDE_DIR
      AP_PPLITE_LIB
      FLINT_FOUND
      GMPXX_FOUND
      GMP_FOUND
    FAIL_MESSAGE
      "Could NOT find PPLite. Please provide -DPPLITE_ROOT=/path/to/pplite")
endif()

set(PPLITE_INCLUDE_DIRS
  ${AP_PPLITE_INCLUDE_DIR}
  ${PPLITE_INCLUDE_DIR}
  ${FLINT_INCLUDE_DIR}
  ${GMPXX_INCLUDE_DIR}
  ${GMP_INCLUDE_DIR})

set(PPLITE_LIBRARY
  ${AP_PPLITE_LIB}
  ${PPLITE_LIB}
  ${FLINT_LIB}
  ${GMPXX_LIB}
  ${GMP_LIB})
