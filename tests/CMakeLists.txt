set (CRAB_LIBS 
  dbm
  term)

if (USE_LDD)
  set (CRAB_LIBS ${CRAB_LIBS} boxes) 
endif ()

set (CRAB_LIBS ${CRAB_LIBS}
  ${LDD_LIBRARY} 
  ${GMPXX_LIB} 
  ${GMP_LIB} 
  ${RT_LIB}
  )

add_subdirectory (simple)
if (USE_LDD)
  add_subdirectory (boxes)
endif ()
add_subdirectory (arrays)
add_subdirectory (pointers)
add_subdirectory (unif)
add_subdirectory (inter)
add_subdirectory (conc)
# add_subdirectory (cg)

