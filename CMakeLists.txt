project(ikos-core)

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

find_package(GMP REQUIRED)
if (GMP_FOUND)
  include_directories (${GMP_INCLUDE_DIR})
  include_directories (${GMPXX_INCLUDE_DIR})
else()
  set(GMP_LIB "")
  set(GMPXX_LIB "")
endif()

include_directories (${CMAKE_SOURCE_DIR})

# check for rt lib. Not needed on OSX.
find_library(RT_LIB NAMES rt)
if (NOT RT_LIB)
  set(RT_LIB "")
endif()
mark_as_advanced(RT_LIB)

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#add_definitions(-std=c++11)
add_definitions(-Wno-redeclared-class-member -Wno-sometimes-uninitialized)
add_definitions(-Wno-covered-switch-default)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
add_definitions( -Wno-unused-local-typedefs)
endif ()

add_subdirectory(ikos/dbm)
add_subdirectory(tests)



