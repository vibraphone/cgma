
PROJECT(exo)

ADD_EXECUTABLE(exo ExoLoader.cpp ElemConnectivity.cpp)

INCLUDE_DIRECTORIES(${CUBIT_SOURCE_BASE}/util ${CUBIT_SOURCE_BASE}/geom
                    ${CUBIT_BUILD_BASE}/util ${CUBIT_BUILD_BASE}/geom
                    ${CUBIT_SOURCE_BASE}/geom/facet
                    ${CUBIT_SOURCE_BASE}/geom/Cholla
                    ${EXODUS_INCLUDE_DIR}
                    ${NETCDF_INCLUDE_DIR})

LINK_DIRECTORIES(${CUBIT_LIB_BASE} ${EXODUS_LIBRARY_DIR}
                 ${NETCDF_LIBRARY_DIR})

TARGET_LINK_LIBRARIES(exo cubit_facet cubit_facetboolstub cholla 
                      cubit_geom cubit_util ${EXODUS_LIBRARY} netcdf)

##If doing an out of source build copy the test
##input file.
IF(NOT ${exo_SOURCE_DIR} STREQUAL ${exo_BINARY_DIR})
  ADD_CUSTOM_COMMAND(OUTPUT ${exo_BINARY_DIR}/exo.g
                     COMMAND ${CMAKE_COMMAND}
                     ARGS -E copy ${exo_SOURCE_DIR}/exo.g ${exo_BINARY_DIR}/exo.g
                     MAIN_DEPENDENCY ${exo_SOURCE_DIR}/exo.g)
  ADD_CUSTOM_TARGET(copy_exo ALL
                    DEPENDS ${exo_BINARY_DIR}/exo.g)
ENDIF(NOT ${exo_SOURCE_DIR} STREQUAL ${exo_BINARY_DIR})

