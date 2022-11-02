# Find FLINT headers and libraries.
# Author: Enea Zaffanella (enea.zaffanella@unipr.it)

if (NOT FLINT_FOUND)
  set(FLINT_ROOT "" CACHE PATH "Path to flint install directory")

  find_path(FLINT_INCLUDES
    NAMES flint/flint.h
    HINTS "${FLINT_ROOT}/include"
    DOC "Path to flint include directory"
  )

  find_library(FLINT_LIB
    NAMES flint
    HINTS "${FLINT_ROOT}/lib"
    DOC "Path to flint library"
  )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(FLINT
    REQUIRED_VARS
      FLINT_INCLUDES
      FLINT_LIB
    VERSION_VAR
      FLINT_VERSION
    FAIL_MESSAGE
      "Could NOT find FLINT. Please provide -DFLINT_ROOT=/path/to/flint")
endif()
