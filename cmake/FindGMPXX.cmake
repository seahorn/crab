if (NOT GMPXX_FOUND)
  set(GMP_SEARCH_PATH "" CACHE PATH "Search path for gmp.")
  find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h PATHS ${GMP_SEARCH_PATH}/include)
  find_library(GMPXX_LIB NAMES gmpxx PATHS ${GMP_SEARCH_PATH}/lib)

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GMPXX
    REQUIRED_VARS GMPXX_INCLUDE_DIR GMPXX_LIB)

endif()
