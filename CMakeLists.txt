project(Crab)

cmake_minimum_required (VERSION 2.8.11)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR ) 
  message (FATAL_ERROR 
    "In-source builds are not allowed. Please clean your source tree and try again.")  
endif()

# Add path for custom modules
set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall")

set (CUSTOM_BOOST_ROOT "" CACHE PATH "Path to custom boost installation.")
if (CUSTOM_BOOST_ROOT)
  set (BOOST_ROOT ${CUSTOM_BOOST_ROOT})
  set (Boost_NO_SYSTEM_PATHS "ON")
endif()

set (Boost_USE_STATIC_LIBS ON)
find_package (Boost REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})

option (DEVMODE   "Enable development mode" OFF)
option (USE_LDD   "Enable Ldd library" OFF)
option (USE_APRON "Enable Apron library" OFF)

find_package(GMP REQUIRED)
if (GMP_FOUND)
  include_directories (${GMP_INCLUDE_DIR})
  include_directories (${GMPXX_INCLUDE_DIR})
else()
  set(GMP_LIB "")
  set(GMPXX_LIB "")
endif()

include(ExternalProject)
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS configure build test)

find_package (Git)

if (USE_LDD)
if (GIT_FOUND)
    set (LDD_TAG "origin/master" CACHE STRING "ldd tag to use")
    set (LDD_REPO "http://git.code.sf.net/p/lindd/code" CACHE STRING "ldd repository")
    ExternalProject_Add(ldd
      GIT_REPOSITORY ${LDD_REPO}
      GIT_TAG ${LDD_TAG}
      PREFIX ${CMAKE_BINARY_DIR}/ldd
      INSTALL_DIR ${CMAKE_BINARY_DIR}/run/ldd
      CMAKE_ARGS
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      TEST_AFTER_INSTALL 1
      TEST_COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_LIST_FILE}
      LOG_DOWNLOAD 1
      LOG_CONFIGURE 1
      LOG_BUILD 1
      LOG_INSTALL 1)
else ()
  message (STATUS "Could not find git. Not downloading ldd")
endif()
endif ()


if (USE_APRON)
if (GIT_FOUND)
  get_filename_component (GMP_SEARCH_PATH ${GMP_INCLUDE_DIR} PATH)

  find_package(MPFR REQUIRED)
  if (MPFR_FOUND)
    get_filename_component (MPFR_SEARCH_PATH ${MPFR_INCLUDES} PATH)
  endif ()

  set (APRON_TAG "origin/elina" CACHE STRING "apron tag to use")
  ExternalProject_Add(apron
    GIT_REPOSITORY git@bitbucket.org:jorgenavas/apron.git 
    GIT_TAG ${APRON_TAG}
    CONFIGURE_COMMAND echo "Apron does not need a configure"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} 
    APRON_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} 
    INSTALL_COMMAND 
    make APRON_PREFIX=<INSTALL_DIR> GMP_PREFIX=${GMP_SEARCH_PATH} MPFR_PREFIX=${MPFR_SEARCH_PATH} install
    LOG_CONFIGURE 1
    LOG_INSTALL 1
    LOG_BUILD 1)
else ()
  message (STATUS "Could not find git. Not downloading apron")
endif()
endif()

if (USE_LDD)
  find_package(Ldd)
  if (LDD_FOUND)
    include_directories (${LDD_INCLUDE_DIR})
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LDD_CXXFLAGS}")
    set (LDD_LIBS ${LDD_LIBRARY})
    set (HAVE_LDD TRUE)
  else()
     ExternalProject_Get_Property (ldd INSTALL_DIR)
     set (LDD_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of ldd" FORCE)
     message (WARNING "No ldd found. Run \n\tcmake --build . --target ldd && cmake ${CMAKE_SOURCE_DIR}")
     return()
  endif()
endif ()

if (USE_APRON)
  find_package(Apron)
  if (APRON_FOUND)
    include_directories (${APRON_INCLUDE_DIR})
    set (APRON_LIBS ${APRON_LIBRARY})
    set (HAVE_APRON TRUE)
  else()
    ExternalProject_Get_Property (apron INSTALL_DIR)
    set (APRON_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of apron" FORCE)
    message (WARNING "No apron found. Run \n\tcmake --build . --target apron && cmake ${CMAKE_SOURCE_DIR}")
    return()
  endif()
endif ()

include_directories (${Crab_SOURCE_DIR}/include)

if (USE_LDD)
  set_target_properties(ldd PROPERTIES EXCLUDE_FROM_ALL ON)
endif ()

if (USE_APRON)
  set_target_properties(apron PROPERTIES EXCLUDE_FROM_ALL ON)
endif ()


# check for rt lib. Not needed on OSX.
find_library(RT_LIB NAMES rt)
if (NOT RT_LIB)
  set(RT_LIB "")
endif()
mark_as_advanced(RT_LIB)

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_definitions(-std=c++11)
add_definitions(-Wno-redeclared-class-member -Wno-sometimes-uninitialized)
add_definitions(-Wno-covered-switch-default)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
add_definitions( -Wno-unused-local-typedefs)
endif ()

include_directories (${CMAKE_BINARY_DIR}/include)

configure_file( include/crab/config.h.cmake
                ${CMAKE_BINARY_DIR}/include/crab/config.h )

add_subdirectory(lib)

if (DEVMODE)
   message (STATUS "Enabled development mode: tests will be compiled")
   add_subdirectory(tests)
endif ()

install(DIRECTORY include/ DESTINATION include)

# Licenses
install(FILES Crab_LICENSE.txt DESTINATION shared/doc)       

install(FILES Ikos_LICENSE.pdf DESTINATION shared/doc)       





