
PROJECT(facets)


ADD_EXECUTABLE(facets facets.cpp)

INCLUDE_DIRECTORIES(${CUBIT_SOURCE_BASE}/util ${CUBIT_SOURCE_BASE}/geom
                    ${CUBIT_BUILD_BASE}/util ${CUBIT_BUILD_BASE}/geom
                    ${CUBIT_SOURCE_BASE}/geom/facet
                    ${CUBIT_SOURCE_BASE}/geom/Cholla )

LINK_DIRECTORIES(${CUBIT_LIB_BASE})

TARGET_LINK_LIBRARIES(facets cubit_facet cubit_facetboolstub cholla cubit_geom cubit_util)

ADD_CUSTOM_COMMAND(TARGET facets
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND}
		   ARGS -E copy ${facets_BINARY_DIR}/facets ${facets_BINARY_DIR}/facets2)
