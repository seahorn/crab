#Find Elina library
if (NOT ELINA_FOUND)

   set (ELINA_ROOT "" CACHE PATH "Root of Elina install.")

   find_package (Gmp QUIET)
   find_package (Mpfr QUIET)

   find_path(ELINA_INCLUDE_DIR NAMES elina_abstract0.h PATHS ${ELINA_ROOT}/include)

   ### XXX: we grab static libraries
   # find_library(Elina_Poly_Lib NAMES liboptpoly.a PATHS ${ELINA_ROOT}/lib)
   # find_library(Elina_Opt_Oct_Lib NAMES liboptoct.a PATHS ${ELINA_ROOT}/lib)
   # find_library(Elina_Linearize_Lib NAMES libelinalinearize.a PATHS ${ELINA_ROOT}/lib)
   # find_library(Elina_Aux_Lib NAMES libelinaux.a PATHS ${ELINA_ROOT}/lib)
   # find_library(Elina_Opt_Zones_Lib NAMES liboptzones.a PATHS ${ELINA_ROOT}/lib)
   # find_library(Elina_Zonotope_Lib NAMES libzonotope.a PATHS ${ELINA_ROOT}/lib)      

   ### XXX: we grab dynamic libraries
   find_library(Elina_Poly_Lib NAMES liboptpoly.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Opt_Oct_Lib NAMES liboptoct.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Linearize_Lib NAMES libelinalinearize.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Aux_Lib NAMES libelinaux.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Opt_Zones_Lib NAMES liboptzones.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Zonotope_Lib NAMES libzonotope.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   find_library(Elina_Partitions_Lib NAMES libpartitions.so PATHS ${ELINA_ROOT}/lib NO_DEFAULT_PATH)
   
   set(ELINA_LIBRARY ${Elina_Poly_Lib} 
     ${Elina_Opt_Oct_Lib} ${Elina_Linearize_Lib} ${Elina_Aux_Lib}
     ${Elina_Opt_Zones_Lib} ${Elina_Zonotope_Lib} ${Elina_Partitions_Lib})

   include (FindPackageHandleStandardArgs)
   find_package_handle_standard_args (Elina
     REQUIRED_VARS ELINA_INCLUDE_DIR ELINA_LIBRARY GMP_FOUND MPFR_FOUND)
   
   set (ELINA_INCLUDE_DIR ${ELINA_INCLUDE_DIR} ${MPFR_INC_DIR})
   set (ELINA_LIBRARY ${ELINA_LIBRARY} ${MPFR_LIB})
   
   mark_as_advanced(ELINA_LIBRARY ELINA_INCLUDE_DIR)
endif ()
