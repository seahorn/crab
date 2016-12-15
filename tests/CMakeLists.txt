set (CRAB_LIBS 
   OS
   Term
   Debug
   Stats
   )

set (CRAB_LIBS ${CRAB_LIBS}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${LDD_LIBRARY} 
  ${APRON_LIBRARY}
  ${MPFR_LIBRARIES} 
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
add_subdirectory (nullity)
add_subdirectory (terms)
add_subdirectory (inter)
add_subdirectory (cg)
add_subdirectory (cfg)
add_subdirectory (thresholds)
add_subdirectory (disintervals)
add_subdirectory (checkers)

# Needed if shared libraries
if (BUILD_CRAB_LIBS_SHARED)
if (USE_APRON)
  install(FILES ${APRON_LIBRARY} DESTINATION lib)
endif ()
if (USE_LDD)
  install(FILES ${LDD_LIBRARY} DESTINATION lib)
endif ()
endif ()
