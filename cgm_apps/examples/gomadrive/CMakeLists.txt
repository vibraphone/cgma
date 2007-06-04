
PROJECT(drive)

ADD_EXECUTABLE(drive main.cpp gm_cgm_c_interface.C)

INCLUDE_DIRECTORIES(${drive_SOURCE_DIR} ${CUBIT_SOURCE_BASE}/util 
                    ${CUBIT_SOURCE_BASE}/geom ${CUBIT_BUILD_BASE}/util 
                    ${CUBIT_BUILD_BASE}/geom
                    ${CUBIT_SOURCE_BASE}/geom/ACIS ${ACIS_INCLUDE_DIR}
                    ${CUBIT_SOURCE_BASE}/geom/virtual)

LINK_DIRECTORIES(${CUBIT_LIB_BASE} ${ACIS_LIBRARY_DIR})

TARGET_LINK_LIBRARIES(drive ${ACIS_LIBS} cubit_geom cubit_virtual cubit_facet cubit_facetboolstub cubit_geom cubit_util cubit_ACIS cholla)

ADD_DEFINITIONS( -DUSE_CGM )
