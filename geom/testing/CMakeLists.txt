
# test suite to test cgm apis
# this test suite attempts to work for all cgm ports
PROJECT(CGM_TESTS)


# a list of tests
# each testname has a testname.cpp file that implements 
# a "int testname(int, char**)" as its entry point
# testname(int, char**) returns 0 on success
# when testname(int, char**) is called, cgm has already been initialized
# and is ready to use with active ports already set
SET(CGM_TESTS
    AngleCalc
    CreateGeometry
    ReadIgesFile
    GraphicsData
   )

 
# CGM_PORTS contains an identifier for the port to use
# these identifiers are recognized in cgm_test.cpp
IF(CGM_ACIS)
  SET(CGM_PORTS ${CGM_PORTS} "acis")
  INCLUDE(${cubit_acis_SOURCE_DIR}/UseACIS.cmake)
  ADD_DEFINITIONS(-DCGM_ACIS)
  INCLUDE_DIRECTORIES(${cubit_acis_SOURCE_DIR} ${ACIS_INCLUDE_DIR})
  SET(PORT_LIBS ${PORT_LIBS} cubit_ACIS gtcAttrib)
  IF(NOT WIN32 AND BUILD_64)
    SET(SKIP_ReadIgesFile-acis 1)
  ENDIF(NOT WIN32 AND BUILD_64)
ENDIF(CGM_ACIS)
IF(CGM_SMLIB)
  SET(CGM_PORTS ${CGM_PORTS} "smlib")
  INCLUDE(${cubit_smlib_SOURCE_DIR}/UseSMLib.cmake)
  INCLUDE_DIRECTORIES(${cubit_smlib_SOURCE_DIR} ${SMLIB_INCLUDE_DIR})
  ADD_DEFINITIONS(-DCGM_SMLIB)
  SET(PORT_LIBS ${PORT_LIBS} cubit_smlib)
ENDIF(CGM_SMLIB)


# functions that'll start/end cgm
# CREATE_TEST_SOURCELIST will insert the function calls into the test driver
SET(CMAKE_TESTDRIVER_BEFORE_TESTMAIN "  start_cgm(ac,av);")
SET(CMAKE_TESTDRIVER_AFTER_TESTMAIN "  end_cgm(ac,av);")
CREATE_TEST_SOURCELIST(TEST_SRCS cgm_test_driver.cpp ${CGM_TESTS} 
                       EXTRA_INCLUDE "cgm_test.hpp")

SET(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/TestConfig.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/TestConfig.h @ONLY)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(TEST_SRCS
  ${TEST_SRCS}
  cgm_test.cpp
  cgm_test.hpp
  TestUtilities.cpp
  TestUtilities.hpp
  TestConfig.h
  )

ADD_EXECUTABLE(cgm_test ${TEST_SRCS})
TARGET_LINK_LIBRARIES(cgm_test ${PORT_LIBS} cubit_geom cubit_util)

GET_TARGET_PROPERTY(cgm_test_exe cgm_test LOCATION)

# add all the tests for each port
FOREACH(port ${CGM_PORTS})
  MESSAGE("adding tests for ${port}")
  FOREACH(test ${CGM_TESTS})
    IF(NOT SKIP_${test}-${port})
      ADD_TEST(
        cgm-${test}-${port} ${cgm_test_exe} ${test} 
        -E ${port}
      )
    ENDIF(NOT SKIP_${test}-${port})
  ENDFOREACH(test)
ENDFOREACH(port)

