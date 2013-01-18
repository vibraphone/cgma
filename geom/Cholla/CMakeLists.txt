
PROJECT(Cholla)

SET(CHOLLA_SRCS
    AllocMemManagersCholla.cpp
    Cholla.cpp
    ChollaCurve.cpp
    ChollaCurve.hpp
    ChollaEngine.cpp
    ChollaEngine.hpp
    ChollaEntity.cpp
    ChollaEntity.hpp
    ChollaPoint.cpp
    ChollaPoint.hpp
    ChollaSkinTool.cpp
    ChollaSkinTool.hpp
    ChollaSurface.cpp
    ChollaSurface.hpp
    ChollaVolume.cpp
    ChollaVolume.hpp
    ChordalAxis.cpp
    ChordalAxis.hpp
    CubitFacet.cpp
    CubitFacet.hpp
    CubitFacetData.cpp
    CubitFacetData.hpp
    CubitFacetEdge.cpp
    CubitFacetEdge.hpp
    CubitFacetEdgeData.cpp
    CubitFacetEdgeData.hpp
    CubitPoint.cpp
    CubitPoint.hpp
    CubitPointData.cpp
    CubitPointData.hpp
    CubitQuadFacet.cpp
    CubitQuadFacet.hpp
    CubitQuadFacetData.cpp
    CubitQuadFacetData.hpp
    CurveFacetEvalTool.cpp
    CurveFacetEvalTool.hpp
    debug.cpp
    debug.hpp
    FacetDataUtil.cpp
    FacetDataUtil.hpp
    FacetEntity.cpp
    FacetEntity.hpp
    FacetEvalTool.cpp
    FacetEvalTool.hpp
    GeoNode.cpp
    GeoNode.hpp
    GeoTet.cpp
    GeoTet.hpp
    LoopParamTool.cpp
    LoopParamTool.hpp
    PointGridSearch.cpp
    PointGridSearch.hpp
    PointLoopFacetor.cpp
    PointLoopFacetor.hpp
    TDChordal.cpp
    TDChordal.hpp
    TDFacetboolData.cpp
    TDFacetboolData.hpp
    TDFacetBoundaryEdge.cpp
    TDFacetBoundaryEdge.hpp
    TDFacetBoundaryPoint.cpp
    TDFacetBoundaryPoint.hpp
    TDGeomFacet.cpp
    TDGeomFacet.hpp
    )

#INCLUDE_DIRECTORIES(${Cholla_SOURCE_DIR})
INCLUDE_DIRECTORIES(${cubit_util_SOURCE_DIR})
INCLUDE_DIRECTORIES(${cubit_util_BINARY_DIR})

ADD_LIBRARY(cholla STATIC ${CHOLLA_SRCS})


