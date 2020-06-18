project(Crab)

cmake_minimum_required(VERSION 3.3)

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
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE)
endif()

# Add path for custom modules
set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

#Setup CMake to run tests
enable_testing()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set (CUSTOM_BOOST_ROOT "" CACHE PATH "Path to custom boost installation.")
if (CUSTOM_BOOST_ROOT)
  set (BOOST_ROOT ${CUSTOM_BOOST_ROOT})
  set (Boost_NO_SYSTEM_PATHS "ON")
endif()


option (CRAB_USE_LDD   "Enable Ldd library" OFF)
option (CRAB_USE_APRON "Enable Apron library" OFF)
option (CRAB_USE_ELINA "Enable Elina library" OFF)
option (CRAB_BUILD_LIBS_SHARED "Build all Crab libraries dynamically." OFF)
option (CRAB_ENABLE_TESTS "Enable tests" OFF)
option (CRAB_ENABLE_STATS "Enable stats" ON)

if (CRAB_USE_APRON AND CRAB_USE_ELINA)
  message(FATAL_ERROR "Apron and Elina are not compatible: choose one")
endif()

# Directory where all (binary) tests will be located
# Binary tests do not go to install directory by default
if (CRAB_ENABLE_TESTS)
  if (NOT TEST_DIR)
    set(TEST_DIR ${CMAKE_BINARY_DIR}/test-bin)
  endif()  
endif()

# Save old CMAKE_FIND_LIBRARY_SUFFIXES
set(_CRAB_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
#set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib .a)
if (NOT CRAB_BUILD_LIBS_SHARED)
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif ()   

# To remember the last configuration
if (CRAB_ENABLE_TESTS)
  message(STATUS "cmake crab configuration: cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCRAB_USE_APRON=${CRAB_USE_APRON} -DCRAB_USE_ELINA=${CRAB_USE_ELINA} -DCRAB_USE_LDD=${CRAB_USE_LDD} -DCRAB_BUILD_LIBS_SHARED=${CRAB_BUILD_LIBS_SHARED} -DCRAB_ENABLE_TESTS=${CRAB_ENABLE_TESTS} -DTEST_DIR=${TEST_DIR}")
else()
  message(STATUS "cmake crab configuration: cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCRAB_USE_APRON=${CRAB_USE_APRON} -DCRAB_USE_ELINA=${CRAB_USE_ELINA} -DCRAB_USE_LDD=${CRAB_USE_LDD} -DCRAB_BUILD_LIBS_SHARED=${CRAB_BUILD_LIBS_SHARED}")
endif()  

if (CRAB_BUILD_LIBS_SHARED)
  set(_DYNAMIC_LIBS_EXC_STR "")
  if(CRAB_USE_LDD)
    set(_DYNAMIC_LIBS_EXC_STR "${_DYNAMIC_LIBS_EXC_STR}\n  - except ldd ")
  endif()
  if(CRAB_USE_APRON)
    set(_DYNAMIC_LIBS_EXC_STR "${_DYNAMIC_LIBS_EXC_STR}\n  - except apron")
  endif()
  message (STATUS "Crab libraries and dependencies are built dynamically ${_DYNAMIC_LIBS_EXC_STR}")  
  set(CRAB_LIBS_TYPE SHARED)
else()
  set(_STATIC_LIBS_EXC_STR "")
  if(CRAB_USE_ELINA)
    set(_STATIC_LIBS_EXC_STR "${_STATIC_LIBS_EXC_STR}\n  - except elina")
  endif()
  message (STATUS "Crab libraries and dependencies are built statically ${_STATIC_LIBS_EXC_STR}")  
  set(CRAB_LIBS_TYPE STATIC)
endif ()

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

if (CRAB_ENABLE_TESTS)
  find_package (Boost 1.65.0 COMPONENTS program_options REQUIRED)
  set(Boost_USE_STATIC_LIBS ON)
  if (CRAB_BUILD_LIBS_SHARED)
    set(Boost_USE_STATIC_LIBS OFF)
  endif ()
  add_definitions(-DBOOST_NO_EXCEPTIONS)
else ()
  find_package (Boost 1.65.0)
endif ()
if (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif ()

set(GMP_USE_STATIC_LIBS ON)
if (CRAB_BUILD_LIBS_SHARED)
  set(GMP_USE_STATIC_LIBS OFF)
endif ()
  
find_package(GMP REQUIRED)
if (GMP_FOUND)
  include_directories (${GMP_INCLUDE_DIR})
else()
  set(GMP_LIB "")
endif()

include(ExternalProject)
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS configure build test)

# Configuration for Coverage
if (CMAKE_BUILD_TYPE STREQUAL Coverage)
  include(Coverage)
  # We don't want coverage on external things
  set(EXT_CMAKE_BUILD_TYPE Release)
else()
  set(EXT_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif ()

find_package (Git)

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
  message (STATUS "Could not find git. Not downloading ldd")
endif()
endif ()


if (CRAB_USE_APRON OR CRAB_USE_ELINA)
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
  message (STATUS "Could not find git. Not downloading apron")
endif()
endif()


if (CRAB_USE_ELINA)
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

if (CRAB_USE_LDD)
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
     # restore old CMAKE_FIND_LIBRARY_SUFFIXES
     set(CMAKE_FIND_LIBRARY_SUFFIXES ${_CRAB_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})     
     return()
  endif()
endif ()

if (CRAB_USE_APRON)
  find_package(Apron)
  if (APRON_FOUND)
    include_directories (${APRON_INCLUDE_DIR})
    set (APRON_LIBS ${APRON_LIBRARY})
    set (HAVE_APRON TRUE)
  else()
    ExternalProject_Get_Property (apron INSTALL_DIR)
    set (APRON_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of apron" FORCE)
    message (WARNING "No apron found. Run \n\tcmake --build . --target apron && cmake ${CMAKE_SOURCE_DIR}")
    # restore old CMAKE_FIND_LIBRARY_SUFFIXES
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_CRAB_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})     
    return()
  endif()
endif ()

if (CRAB_USE_ELINA)
  find_package(Elina)
  if (ELINA_FOUND)
    include_directories (${ELINA_INCLUDE_DIR})
    set (ELINA_LIBS ${ELINA_LIBRARY})
    set (HAVE_ELINA TRUE)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")    
      set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags")
    endif()
  else()
    ExternalProject_Get_Property (elina INSTALL_DIR)
    set (ELINA_ROOT ${INSTALL_DIR} CACHE FILEPATH "Forced location of elina" FORCE)
    message (WARNING "No elina found. Run \n\tcmake --build . --target elina && cmake ${CMAKE_SOURCE_DIR}")
    # restore old CMAKE_FIND_LIBRARY_SUFFIXES
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_CRAB_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})         
    return()
  endif()
endif ()

include_directories (${Crab_SOURCE_DIR}/include)

if (CRAB_USE_LDD)
  set_target_properties(ldd PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The boxes domain will be available")  
else ()
  message (STATUS "The boxes domain will not be available")  
endif ()

if (CRAB_USE_APRON)
  set_target_properties(apron PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The apron domains will be available")
else ()
  message (STATUS "The apron domains will not be available")
endif ()

if (CRAB_USE_ELINA)
  set_target_properties(elina PROPERTIES EXCLUDE_FROM_ALL ON)
  message (STATUS "The elina domains will be available")
else ()
  message (STATUS "The elina domains will not be available")
endif ()

# Enable always logging
set(NCRABLOG FALSE)

# Enable or disable stats gathering
if (CRAB_ENABLE_STATS)
  message (STATUS "Compiled Crab with support for gathering statistics")  
  set(CRAB_STATS TRUE)
else()
  message (STATUS "Compiled Crab without support for gathering statistics")    
  set(CRAB_STATS FALSE)
endif()
    
# check for rt lib. Not needed on OSX.
# find_library(RT_LIB NAMES rt)
# if (NOT RT_LIB)
#   set(RT_LIB "")
# endif()
# mark_as_advanced(RT_LIB)

# XXX: needed if shared libraries
# if (CRAB_BUILD_LIBS_SHARED)
#   if (CRAB_USE_APRON)
#     if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
#       install(FILES ${APRON_LIBRARY} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#     else ()
#       install(FILES ${APRON_LIBRARY} DESTINATION lib)
#     endif ()
#   endif ()
#   if (CRAB_USE_LDD)
#     if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
#       install(FILES ${LDD_LIBRARY} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#     else ()
#       install(FILES ${LDD_LIBRARY} DESTINATION lib)
#     endif ()
#   endif ()
# endif ()


# XXX: do not use add_definition for warning flags. add_definition is
#      intended for -D flags.
###
# To avoid annoying warning messages
###
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
else ()
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-redeclared-class-member")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-covered-switch-default")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sometimes-uninitialized")    
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-uninitialized")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef")
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-inconsistent-missing-override")
  endif ()
endif ()

set (CRAB_DEPS_LIBS
  ${LDD_LIBRARY} 
  ${APRON_LIBRARY}
  ${ELINA_LIBRARY}
  ${MPFR_LIBRARIES} 
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
  if (ELINA_FOUND)
     set (_CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} ${ELINA_INCLUDE_DIR})
  endif ()
   
  set (CRAB_INCLUDE_DIRS ${_CRAB_INCLUDE_DIRS} PARENT_SCOPE)  
endif ()  
  
include_directories (${CMAKE_BINARY_DIR}/include)

configure_file(include/crab/config.h.cmake
                ${CMAKE_BINARY_DIR}/include/crab/config.h )

add_subdirectory(lib)

if (CRAB_ENABLE_TESTS)
   message (STATUS "Tests will be compiled")
   add_subdirectory(tests)
endif ()

# restore old CMAKE_FIND_LIBRARY_SUFFIXES
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_CRAB_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})     

install(DIRECTORY include/
        DESTINATION crab/include
        PATTERN "config.h.cmake" EXCLUDE)
install(FILES ${CMAKE_BINARY_DIR}/include/crab/config.h DESTINATION crab/include/crab)

if (CRAB_USE_LDD)
  ExternalProject_Get_Property (ldd INSTALL_DIR)  
  install(DIRECTORY ${INSTALL_DIR}
         DESTINATION .)
endif()

if (CRAB_USE_APRON)
  ExternalProject_Get_Property (apron INSTALL_DIR)  
  install(DIRECTORY ${INSTALL_DIR}
          DESTINATION .)
endif()

if (CRAB_USE_ELINA)
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






