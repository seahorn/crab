set (CRAB_LIBS 
  dbm
  term
  ${GMPXX_LIB}
  ${GMP_LIB}
  )

add_subdirectory (simple)
add_subdirectory (arrays)
add_subdirectory (pointers)
add_subdirectory (unif)
add_subdirectory (conc)
