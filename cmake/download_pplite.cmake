if (CRAB_USE_PPLITE)
  ExternalProject_Add(pplite
    URL https://github.com/ezaffanella/PPLite/raw/main/releases/pplite-0.9.tar.gz
    INSTALL_DIR ${CMAKE_BINARY_DIR}/run/pplite
    CONFIGURE_COMMAND
    ../pplite/configure --prefix <INSTALL_DIR> --enable-apron --with-gmp=${GMP_SEARCH_PATH} --with-flint=${FLINT_SEARCH_PATH} --with-apron=${APRON_ROOT}
    BUILD_IN_SOURCE 0
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    LOG_CONFIGURE 1
    LOG_INSTALL 1
    LOG_BUILD 1)
endif()
