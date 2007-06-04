
PROJECT(mergechk)


ADD_EXECUTABLE(mergechk mergechk.cpp)

INCLUDE_DIRECTORIES(${CUBIT_SOURCE_BASE}/util ${CUBIT_SOURCE_BASE}/geom
                    ${CUBIT_BUILD_BASE}/util ${CUBIT_BUILD_BASE}/geom
                    ${CUBIT_SOURCE_BASE}/geom/ACIS 
                    ${DUBIT_SOURCE_BASE}/geom/Cholla ${ACIS_INCLUDE_DIR})

LINK_DIRECTORIES(${CUBIT_LIB_BASE} ${ACIS_LIBRARY_DIR})

TARGET_LINK_LIBRARIES(mergechk ${ACIS_LIBS} cubit_geom cubit_facet cubit_facetboolstub cubit_geom cubit_util cubit_ACIS cholla cubit_geom)

##If doing an out of source build copy the test
##input file.
IF(NOT ${mergechk_SOURCE_DIR} STREQUAL ${mergechk_BINARY_DIR})
  ADD_CUSTOM_COMMAND(OUTPUT ${mergechk_BINARY_DIR}/merge.sat
                     COMMAND ${CMAKE_COMMAND}
                     ARGS -E copy ${mergechk_SOURCE_DIR}/merge.sat ${mergechk_BINARY_DIR}/merge.sat
                     MAIN_DEPENDENCY ${mergechk_SOURCE_DIR}/merge.sat)
  ADD_CUSTOM_TARGET(copy_merge ALL
                    DEPENDS ${mergechk_BINARY_DIR}/merge.sat)
ENDIF(NOT ${mergechk_SOURCE_DIR} STREQUAL ${mergechk_BINARY_DIR})
