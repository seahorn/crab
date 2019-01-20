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

#Setup CMake to run tests
enable_testing()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set (CUSTOM_BOOST_ROOT "" CACHE PATH "Path to custom boost installation.")
if (CUSTOM_BOOST_ROOT)
  set (BOOST_ROOT ${CUSTOM_BOOST_ROOT})
  set (Boost_NO_SYSTEM_PATHS "ON")
endif()

option (ENABLE_TESTS "Enable tests" OFF)
option (USE_LDD   "Enable Ldd library" OFF)
option (USE_APRON "Enable Apron library" OFF)
option (USE_MDD   "Enable Mdd library" OFF)
option (USE_ELINA "Enable Elina library" OFF)
option (USE_STATS "Enable statistics" ON)
option (USE_LOG   "Enable logging for debugging purposes" ON)
option (BUILD_CRAB_LIBS_SHARED "Build all Crab libraries dynamically." ON)

if (USE_APRON AND USE_ELINA)
  message(FATAL_ERROR "Apron and Elina are not compatible: choose one")
endif()

# Directory where all (binary) tests will be located
# Binary tests do not go to install directory by default
if (ENABLE_TESTS)
  if (NOT TEST_DIR)
    set(TEST_DIR ${CMAKE_BINARY_DIR}/test-bin)
  endif()  
endif()

if (NOT BUILD_CRAB_LIBS_SHARED)
   # put static libraries first unless we want dynamic libraries
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif ()   

# To remember the last configuration
message(STATUS "cmake configuration: cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DUSE_APRON=${USE_APRON} -DUSE_ELINA=${USE_ELINA} -DUSE_LDD=${USE_LDD} -DUSE_MDD=${USE_MDD} -DENABLE_TESTS=${ENABLE_TESTS} -DTEST_DIR=${TEST_DIR}")

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
if (NOT CMAKE_INSTALL_RPATH )
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/crab/lib")
endif()

   
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
    set (LDD_REPO "https://github.com/seahorn/ldd.git"
         CACHE STRING "ldd repository")
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


if (USE_APRON OR USE_ELINA)
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
endif ()

if (USE_APRON)
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
  message (STATUS "Could not find git. Not downloading apron")
endif()
endif()

if (USE_MDD)
if (GIT_FOUND)
    set (MDD_TAG "origin/master" CACHE STRING "mdd tag to use")
    set (MDD_REPO "https://bitbucket.org/gkgange/mdd-boxes.git"
         CACHE STRING "mdd repository")
    ExternalProject_Add(mdd
      GIT_REPOSITORY ${MDD_REPO}
      GIT_TAG ${MDD_TAG}
      PREFIX ${CMAKE_BINARY_DIR}/mdd
      INSTALL_DIR ${CMAKE_BINARY_DIR}/run/mdd
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/mdd/src/mdd <INSTALL_DIR>
      LOG_DOWNLOAD 1
      LOG_INSTALL 1)
else ()
  message (STATUS "Could not find git. Not downloading mdd")
endif()
endif ()

if (USE_ELINA)
if (GIT_FOUND)
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
  message (STATUS "Could not find git. Not downloading elina")
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

if (USE_MDD)
  find_package(Mdd)
  if (MDD_FOUND)
    include_directories (${MDD_INCLUDE_DIR})
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MDD_CXXFLAGS}")
    set (MDD_LIBS ${MDD_LIBRARY})
    set (HAVE_MDD TRUE)
  else()
     ExternalProject_Get_Property (mdd INSTALL_DIR)
     set (MDD_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of mdd" FORCE)
     message (WARNING "No mdd found. Run \n\tcmake --build . --target mdd && cmake ${CMAKE_SOURCE_DIR}")
     return()
   endif()
 endif()
 
if (USE_ELINA)
  find_package(Elina)
  if (ELINA_FOUND)
    include_directories (${ELINA_INCLUDE_DIR})
    set (ELINA_LIBS ${ELINA_LIBRARY})
    set (HAVE_ELINA TRUE)
  else()
    ExternalProject_Get_Property (elina INSTALL_DIR)
    set (ELINA_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of elina" FORCE)
    message (WARNING "No elina found. Run \n\tcmake --build . --target elina && cmake ${CMAKE_SOURCE_DIR}")
    return()
  endif()
endif ()

include_directories (${Crab_SOURCE_DIR}/include)

if (USE_LDD)
  set_target_properties(ldd PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The ldd-boxes domain will be available")  
else ()
  message (STATUS "The ldd-boxes domain will not be available")  
endif ()

if (USE_APRON)
  set_target_properties(apron PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The apron domains will be available")
else ()
  message (STATUS "The apron domains will not be available")
endif ()

if (USE_MDD)
  set_target_properties(mdd PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The mdd-boxes domain will be available")  
else ()
  message (STATUS "The mdd-boxes domain will not be available")  
endif ()

if (USE_ELINA)
  set_target_properties(elina PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The elina domains will be available")
else ()
  message (STATUS "The elina domains will not be available")
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
# find_library(RT_LIB NAMES rt)
# if (NOT RT_LIB)
#   set(RT_LIB "")
# endif()
# mark_as_advanced(RT_LIB)

# XXX: needed if shared libraries
# if (BUILD_CRAB_LIBS_SHARED)
#   if (USE_APRON)
#     if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
#       install(FILES ${APRON_LIBRARY} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#     else ()
#       install(FILES ${APRON_LIBRARY} DESTINATION lib)
#     endif ()
#   endif ()
#   if (USE_LDD)
#     if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
#       install(FILES ${LDD_LIBRARY} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#     else ()
#       install(FILES ${LDD_LIBRARY} DESTINATION lib)
#     endif ()
#   endif ()
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
  ${MDD_LIBRARY}   
  ${ELINA_LIBRARY}
  ${MPFR_LIBRARIES} 
  ${GMPXX_LIB} 
  ${GMP_LIB})

if (TopLevel)
  set (CRAB_LIBS Crab ${CRAB_DEPS_LIBS})
else ()
  ## XXX: propagate CRAB_LIBS and CRAB_INCLUDE_DIRS to parent
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
  if (MDD_FOUND)
    set (_CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} ${MDD_INCLUDE_DIR})
  endif ()    
  if (ELINA_FOUND)
     set (_CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} ${ELINA_INCLUDE_DIR})
  endif ()
   
  set (CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} PARENT_SCOPE)  
endif ()  
  
include_directories (${CMAKE_BINARY_DIR}/include)

configure_file(include/crab/config.h.cmake
                ${CMAKE_BINARY_DIR}/include/crab/config.h )

add_subdirectory(lib)

if (ENABLE_TESTS)
   message (STATUS "Tests will be compiled")
   add_subdirectory(tests)
endif ()

install(DIRECTORY include/
        DESTINATION crab/include
        PATTERN "config.h.cmake" EXCLUDE)
install(FILES ${CMAKE_BINARY_DIR}/include/crab/config.h DESTINATION crab/include/crab)

if (USE_LDD)
  ExternalProject_Get_Property (ldd INSTALL_DIR)  
  install(DIRECTORY ${INSTALL_DIR}
         DESTINATION .)
endif()

if (USE_APRON)
  ExternalProject_Get_Property (apron INSTALL_DIR)  
  install(DIRECTORY ${INSTALL_DIR}
          DESTINATION .)
endif()

if (USE_MDD)
  ExternalProject_Get_Property (mdd INSTALL_DIR)
  install(DIRECTORY ${INSTALL_DIR}
          DESTINATION .)
endif()
  
if (USE_ELINA)
  ExternalProject_Get_Property (elina INSTALL_DIR)  
  install(DIRECTORY ${INSTALL_DIR}
    DESTINATION .)

  if (NOT TopLevel)
    ## This is ugly fix so that seahorn and crab-llvm can find elina's
    ## shared libraries. It has to do with how CMAKE_INSTALL_RPATH is
    ## handled by these tools.
    foreach(library ${ELINA_LIBRARY})
      get_filename_component(library_name ${library} NAME_WE)
      file(COPY ${library} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endforeach(library)	
  endif()
endif ()

if (TopLevel)
   install(FILES LICENSE DESTINATION .)
   install(FILES README.md DESTINATION .)
   install(FILES IKOS_LICENSE.pdf DESTINATION .)
endif () 






