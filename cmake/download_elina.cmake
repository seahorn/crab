if (CRAB_USE_ELINA)
if (GIT_FOUND)
 find_package(SSE)
  if (AVX_FOUND)
    set (IS_VECTOR "IS_VECTOR=-DVECTOR")
  else ()
    if (SSE2_FOUND OR SSE3_FOUND OR SSSE3_FOUND OR SSE4_1_FOUND)
      set (IS_VECTOR "IS_VECTOR=-DVECTOR")
      set (IS_SSE "IS_SSE=-DSSE")
      message(WARNING "Crab: building elina oct with vectorized operations.")      
    else ()
      message(WARNING "Crab: building elina oct without vectorized operations.")
    endif ()
  endif ()
    
  set (ELINA_TAG "master" CACHE STRING "elina tag to use")
  ExternalProject_Add(elina
    GIT_REPOSITORY https://github.com/eth-sri/ELINA.git
    GIT_TAG ${ELINA_TAG}
    INSTALL_DIR ${CMAKE_BINARY_DIR}/run/elina
    CONFIGURE_COMMAND ./configure -prefix ${INSTALL_DIR} -use-vector  -gmp-prefix ${GMP_SEARCH_PATH} -mpfr-prefix ${MPFR_SEARCH_PATH} 
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${IS_VECTOR} ${IS_SSE}
    ELINA_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} 
    INSTALL_COMMAND 
    make ELINA_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} install
    LOG_CONFIGURE 1
    LOG_INSTALL 1
    LOG_BUILD 1)
else ()
  message (STATUS "Crab: could not find git. Not downloading elina")
endif()
endif()
