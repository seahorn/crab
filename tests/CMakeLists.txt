set (CRAB_LIBS 
  dbm
  term
  ${LDD_LIBRARY} 
  ${GMPXX_LIB} 
  ${GMP_LIB} 
  ${RT_LIB}
  )

add_subdirectory (simple)
add_subdirectory (arrays)
add_subdirectory (pointers)
add_subdirectory (unif)
add_subdirectory (inter)
add_subdirectory (conc)
add_subdirectory (cg)
if (USE_LDD)
  add_subdirectory (boxes)
endif ()

