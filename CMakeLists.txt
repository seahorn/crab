project(Crab)

cmake_minimum_required (VERSION 2.8.11)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR ) 
  message (FATAL_ERROR 
    "In-source builds are not allowed. Please clean your source tree and try again.")  
endif()

# determine if this is top-level or embedded project
if (PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME)
  set (TopLevel TRUE)
else()
  set (TopLevel FALSE)
endif()

# Default is release with debug info
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Add path for custom modules
set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set (CUSTOM_BOOST_ROOT "" CACHE PATH "Path to custom boost installation.")
if (CUSTOM_BOOST_ROOT)
  set (BOOST_ROOT ${CUSTOM_BOOST_ROOT})
  set (Boost_NO_SYSTEM_PATHS "ON")
endif()

option (ENABLE_TESTS "Enable tests" OFF)
option (USE_LDD   "Enable Ldd library" OFF)
option (USE_APRON "Enable Apron library" OFF)
option (USE_STATS "Enable statistics" ON)
option (USE_LOG   "Enable logging for debugging purposes" ON)
option (BUILD_CRAB_LIBS_SHARED "Build all Crab libraries shared." ON)

if (NOT BUILD_CRAB_LIBS_SHARED)
   # put static libraries first unless we want shared libraries
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif ()   

if (BUILD_CRAB_LIBS_SHARED)
  message (STATUS "All crab libraries will be built dynamically")
  set(CRAB_LIBS_TYPE SHARED)
else()
  message (STATUS "All crab libraries will be built statically")
  set(CRAB_LIBS_TYPE STATIC)
endif ()

######
# WARNING: crab dependencies (apron and ldd) are currently compiled
# statically regardless BUILD_CRAB_LIBS_SHARED.
######

# So that executables outside the build tree can find later dynamic
# libraries.  Even if CRAB_LIBS_TYPE=STATIC we could have a mix of
# dynamic and static libraries
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Mac OS X specific code
  set(CMAKE_MACOSX_RPATH TRUE)
endif ()  
if (NOT CMAKE_INSTALL_RPATH)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif ()

if (NOT BUILD_CRAB_LIBS_SHARED)
   set (Boost_USE_STATIC_LIBS ON)
endif ()  
if (ENABLE_TESTS)
  find_package (Boost 1.55.0 COMPONENTS program_options REQUIRED)
else ()
  find_package (Boost 1.55.0)
endif ()
if (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif ()
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
    set (LDD_REPO "https://github.com/seahorn/ldd.git" CACHE STRING "ldd repository")
    ExternalProject_Add(ldd
      GIT_REPOSITORY ${LDD_REPO}
      GIT_TAG ${LDD_TAG}
      PREFIX ${CMAKE_BINARY_DIR}/ldd
      INSTALL_DIR ${CMAKE_BINARY_DIR}/run/ldd
      CMAKE_ARGS
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
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

  find_package(SSE)
  if (AVX_FOUND)
     message(STATUS "Building opt-oct with AVX")
     set (IS_VECTOR "IS_VECTOR=-DVECTOR")
  else ()
     if (SSE2_FOUND OR SSE3_FOUND OR SSSE3_FOUND OR SSE4_1_FOUND)
        message(STATUS "Building opt-oct with SSE")
        set (IS_VECTOR "IS_VECTOR=-DVECTOR")
        set (IS_SSE "IS_SSE=-DSSE")
     else ()
        message(STATUS "Building opt-oct without vectorized operations.")
     endif ()
  endif ()

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
  message (STATUS "The boxes domain will be available")  
else ()
  message (STATUS "The boxes domain will not be available")  
endif ()

if (USE_APRON)
  set_target_properties(apron PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The apron domains will be available")
else ()
  message (STATUS "The apron domains will not be available")
endif ()

if (NOT USE_STATS)
  set (HAVE_STATS FALSE)
  message (STATUS "Analysis statistics is disabled")
else ()
  set (HAVE_STATS TRUE)
  message (STATUS "Analysis statistics is enabled")
endif ()

if (NOT USE_LOG)
  set (NCRABLOG TRUE)
  message (STATUS "Analysis logging for debugging is disabled")
else ()
  set (NCRABLOG FALSE)
  message (STATUS "Analysis logging for debugging is enabled")
endif ()

# check for rt lib. Not needed on OSX.
#find_library(RT_LIB NAMES rt)
#if (NOT RT_LIB)
#  set(RT_LIB "")
#endif()
#mark_as_advanced(RT_LIB)

# find_package(OpenMP)
# if (OpenMP_FOUND)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif ()

# XXX: do not use add_definition for warning flags. add_definition is
#      intended for -D...
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-redeclared-class-member")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-covered-switch-default")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sometimes-uninitialized")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
else ()
   if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-uninitialized")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-inconsistent-missing-override")
   endif ()
endif ()

set (CRAB_DEPS_LIBS
  ${LDD_LIBRARY} 
  ${APRON_LIBRARY}
  ${MPFR_LIBRARIES} 
  ${GMPXX_LIB} 
  ${GMP_LIB})

if (TopLevel)
  set (CRAB_LIBS Crab ${CRAB_DEPS_LIBS})
else ()
  ## XXX: propagate crab libs and header files to parent
  set (CRAB_LIBS Crab ${CRAB_DEPS_LIBS} PARENT_SCOPE)
  set (_CRAB_INCLUDE_DIRS 
       ${Crab_SOURCE_DIR}/include
       ## crab config.h
       ${CMAKE_BINARY_DIR}/include )
  if (LDD_FOUND)
     set (_CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} ${LDD_INCLUDE_DIR})
  endif ()
  if (APRON_FOUND)
     set (_CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} ${APRON_INCLUDE_DIR})
  endif ()
  set (CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} PARENT_SCOPE)
endif ()  

include_directories (${CMAKE_BINARY_DIR}/include)

configure_file( include/crab/config.h.cmake
                ${CMAKE_BINARY_DIR}/include/crab/config.h )

add_subdirectory(lib)

if (ENABLE_TESTS)
   message (STATUS "Tests will be compiled")
   add_subdirectory(tests)
endif ()

install(DIRECTORY include/ DESTINATION include)

if (TopLevel)
   install(FILES LICENSE DESTINATION .)
   install(FILES README.md DESTINATION .)
endif () 





