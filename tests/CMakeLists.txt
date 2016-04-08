set (CRAB_LIBS 
#  dbm
   term
   Debug
   )

set (CRAB_LIBS ${CRAB_LIBS}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${LDD_LIBRARY} 
  ${APRON_LIBRARY}
  ${GMPXX_LIB} 
  ${GMP_LIB} 
  ${RT_LIB}
  )

add_subdirectory (simple)
if (USE_LDD)
  add_subdirectory (boxes)
endif ()
if (USE_APRON)
  add_subdirectory (apron)
endif ()
add_subdirectory (arrays)
add_subdirectory (pointers)
add_subdirectory (nullity)
add_subdirectory (unif)
add_subdirectory (inter)
add_subdirectory (conc)
add_subdirectory (cg)
add_subdirectory (thresholds)
add_subdirectory (disintervals)

