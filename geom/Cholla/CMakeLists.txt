
PROJECT(Cholla)

SET(CHOLLA_SRCS
    ChollaCurve.cpp
    ChollaEngine.cpp
    ChollaEntity.cpp
    ChollaPoint.cpp
    ChollaSkinTool.cpp
    ChollaSurface.cpp
    ChollaVolume.cpp
    ChordalAxis.cpp
    CubitFacet.cpp
    CubitFacetData.cpp
    CubitFacetEdge.cpp
    CubitFacetEdgeData.cpp
    CubitPoint.cpp
    CubitPointData.cpp
    CubitQuadFacet.cpp
    CubitQuadFacetData.cpp
    CurveFacetEvalTool.cpp
    debug.cpp
    FacetDataUtil.cpp
    FacetEntity.cpp
    FacetEvalTool.cpp
    GeoNode.cpp
    GeoTet.cpp
    LoopParamTool.cpp
    PointGridSearch.cpp
    PointLoopFacetor.cpp
    TDChordal.cpp
    TDFacetboolData.cpp
    TDFacetBoundaryEdge.cpp
    TDFacetBoundaryPoint.cpp
    TDGeomFacet.cpp
    )

SET(EXTRA_CHOLLA_SRCS
    AllocMemManagersCholla.cpp
    Cholla.cpp
    )

IF(TEMPLATE_DEFS_INCLUDED)
SET(TEMPLATE_SRCS
    FacetorTool.cpp
    FacetorUtil.cpp
    TetFacetorTool.cpp
    TDDelaunay.cpp
    )
ENDIF(TEMPLATE_DEFS_INCLUDED)
    
IF(WIN32)
  FOREACH(var ${CHOLLA_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(CHOLLA_HEADERS ${CHOLLA_HEADERS} ${header})
  ENDFOREACH(var ${CHOLLA_SRCS})
ENDIF(WIN32)

INCLUDE_DIRECTORIES(${Cholla_SOURCE_DIR})
INCLUDE_DIRECTORIES(${cubit_util_SOURCE_DIR})

ADD_LIBRARY(cholla ${CHOLLA_SRCS} ${EXTRA_CHOLLA_SRCS} ${CHOLLA_HEADERS}
${TEMPLATE_SRCS})


