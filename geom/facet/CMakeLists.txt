
PROJECT(cubit_facet)

SET(FACET_SRCS
    FacetAttrib.cpp
    FacetAttribSet.cpp
    FacetBody.cpp
    FacetboolInterface.cpp
    FacetCoEdge.cpp
    FacetCurve.cpp
    FacetLoop.cpp
    FacetLump.cpp
    FacetModifyEngine.cpp
    FacetParamTool.cpp
    FacetPoint.cpp
    FacetQueryEngine.cpp
    FacetShell.cpp
    FacetSurface.cpp
    GridSearchTree.cpp
    )
    
IF(WIN32)
  FOREACH(var ${FACET_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(FACET_HEADERS ${FACET_HEADERS} ${header})
  ENDFOREACH(var ${FACET_SRCS})
ENDIF(WIN32)

INCLUDE_DIRECTORIES(${cubit_geom_SOURCE_DIR}/facetbool
                    ${cubit_geom_SOURCE_DIR}/Cholla)

ADD_LIBRARY(cubit_facet STATIC ${FACET_SRCS} ${FACET_HEADERS})


