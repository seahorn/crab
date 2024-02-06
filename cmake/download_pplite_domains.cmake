if (CRAB_USE_PPLITE_NATIVE)
  set (PPLITE_DOMAINS_REPO "https://github.com/seahorn/crab-pplite.git"
    CACHE STRING "source of wrapper for pplite domains using native interface")
  set (PPLITE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/crab/domains/pplite")        
  set (PPLITE_WRAPPER "${PPLITE_SOURCE_DIR}/pplite_native_wrapper.hpp")
  
  if (NOT EXISTS ${PPLITE_WRAPPER})
    add_custom_target(pplite-domains-git
      ${GIT_EXECUTABLE} clone ${PPLITE_DOMAINS_REPO} ${PPLITE_SOURCE_DIR})
  endif()
endif()
