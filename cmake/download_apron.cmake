if (CRAB_USE_APRON)
if (GIT_FOUND)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Mac OS X specific code
    set (APRON_TAG "origin/elina-el-capitan" CACHE STRING "apron tag to use")
  else ()
    set (APRON_TAG "origin/elina" CACHE STRING "apron tag to use")
  endif()
  ExternalProject_Add(apron
    GIT_REPOSITORY https://github.com/seahorn/elina.git
    GIT_TAG ${APRON_TAG}
    INSTALL_DIR ${CMAKE_BINARY_DIR}/run/apron
    CONFIGURE_COMMAND echo "Apron does not need a configure"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${IS_VECTOR} ${IS_SSE}
    APRON_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} 
    INSTALL_COMMAND 
    make APRON_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} install
    LOG_CONFIGURE 1
    LOG_INSTALL 1
    LOG_BUILD 1)
else ()
  message (STATUS "Crab: could not find git. Not downloading apron")
endif()
endif()
