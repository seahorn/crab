# Find MDD library
if (NOT MDD_FOUND)

  set(MDD_ROOT "" CACHE PATH "Root of MDD compiled source tree.")

  find_path(MDD_INCLUDE_DIR NAMES include/MDD.hpp PATHS ${MDD_ROOT})
  set(MDD_LIBRARY "")
  
  include (FindPackageHandleStandardArgs)
  # Not required MDD_LIBRARY since for now there is no libraries
  find_package_handle_standard_args(Mdd REQUIRED_VARS MDD_INCLUDE_DIR)
    
  if (MDD_FOUND)  
     set(MDD_CXXFLAGS "-march=native -D__STDC_LIMIT_MACROS")
     set(MDD_LIBRARY ${MDD_LIBRARY})
     set(MDD_INCLUDE_DIR ${MDD_INCLUDE_DIR})
     mark_as_advanced (MDD_INCLUDE_DIR MDD_LIBRARY MDD_CXXFLAGS)
  endif () 

endif ()
