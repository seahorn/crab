add_library(dbm STATIC dbm.cpp expr.cpp)
#target_link_libraries (dbm z)

add_definitions(-Wno-int-to-pointer-cast)
add_definitions(-Wno-sign-compare)

install(TARGETS dbm
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
