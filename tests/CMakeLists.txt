set (IKOS_CORE_LIBS 
  dbm
  term
  ${GMPXX_LIB}
  ${GMP_LIB}
  )

add_subdirectory (conc)
add_subdirectory (simple)
add_subdirectory (arrays)
add_subdirectory (pointers)
add_subdirectory (unif)
