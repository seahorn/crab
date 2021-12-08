if (CRAB_USE_LDD)
if (GIT_FOUND)
    set (LDD_TAG "origin/master" CACHE STRING "ldd tag to use")
    set (LDD_REPO "https://github.com/seahorn/ldd.git"
         CACHE STRING "ldd repository")
    ExternalProject_Add(ldd
      GIT_REPOSITORY ${LDD_REPO}
      GIT_TAG ${LDD_TAG}
      PREFIX ${CMAKE_BINARY_DIR}/ldd
      INSTALL_DIR ${CMAKE_BINARY_DIR}/run/ldd
      CMAKE_ARGS
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE:STRING=${EXT_CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      # XXX: ldd is always compiled statically.
      # We add -fPIC flag so it can be linked with a shared library
      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
      TEST_AFTER_INSTALL 1
      TEST_COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_LIST_FILE}
      LOG_DOWNLOAD 1
      LOG_CONFIGURE 1
      LOG_BUILD 1
      LOG_INSTALL 1)
else ()
  message (STATUS "Crab: could not find git. Not downloading ldd")
endif()
endif ()
