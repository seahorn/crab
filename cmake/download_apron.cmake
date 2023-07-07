if (CRAB_USE_APRON)
  if (GIT_FOUND)
    
    if (CRAB_USE_PPLITE)
      SET(AP_PPLITE_OPTS -pplite-prefix ${PPLITE_ROOT} -flint-prefix ${FLINT_SEARCH_PATH})
    else()
      SET(AP_PPLITE_OPTS -no-pplite)
    endif()

    ExternalProject_Add(apron
      GIT_REPOSITORY https://github.com/antoinemine/apron.git
      GIT_TAG e03832465bdca1888c56ecbe14dcdac0a243dce2
      INSTALL_DIR ${CMAKE_BINARY_DIR}/run/apron
      CONFIGURE_COMMAND 
      ./configure -prefix <INSTALL_DIR> -no-java -no-ocaml -no-ppl ${AP_PPLITE_OPTS} -gmp-prefix ${GMP_SEARCH_PATH} -mpfr-prefix ${MPFR_SEARCH_PATH}
      BUILD_IN_SOURCE 1
      BUILD_COMMAND make    
      INSTALL_COMMAND make install
      LOG_CONFIGURE 1
      LOG_INSTALL 1
      LOG_BUILD 1)
  else ()
    message (STATUS "Crab: could not find git. Not downloading apron")
  endif()
endif()
