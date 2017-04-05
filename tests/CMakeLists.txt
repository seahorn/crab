set (CRAB_LIBS Crab
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${LDD_LIBRARY} 
  ${APRON_LIBRARY}
  ${MPFR_LIBRARIES} 
  ${GMPXX_LIB} 
  ${GMP_LIB} 
  )

function (AddTest testSrc)
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(testName ${testSrc} NAME_WE)
  #Add compile target
  add_executable(${testName} ${testSrc})
  #link with crab libraries
  target_link_libraries(${testName} ${CRAB_LIBS})
  #move testing binaries into a test-bin directory
  set_target_properties(${testName} PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/test-bin)
  #Finally add it to test execution
  add_test(NAME ${testName} 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test-bin
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test-bin/${testName} )
endfunction()

function (AddTestDir dirName)
  file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${dirName}/*.cc)
  foreach(testSrc ${TEST_SRCS})
    AddTest (${testSrc})
  endforeach(testSrc)
endfunction ()  

AddTestDir (domains)
AddTestDir (inter)
AddTestDir (cg)
AddTestDir (cfg)
AddTestDir (thresholds)
AddTestDir (checkers)
