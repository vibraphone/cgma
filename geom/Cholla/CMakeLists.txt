
PROJECT(Cholla)

SET(CHOLLA_SRCS
  AllocMemManagersCholla.cpp
  Cholla.cpp
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

SET(CHOLLA_HEADERS
  BoundaryConstrainTool.cpp
  BoundaryConstrainTool.hpp
  Cholla.h
  ChollaPoint.hpp
  ChollaSkinTool.hpp
  ChollaSurface.hpp
  ChollaVolume.hpp
  ChordalAxis.hpp
  ChollaCurve.hpp
  ChollaEngine.hpp
  ChollaEntity.hpp
  CubitFacet.hpp
  CubitFacetData.hpp
  CubitFacetEdge.hpp
  CubitFacetEdgeData.hpp
  CubitPoint.hpp
  CubitPointData.hpp
  CubitQuadFacet.hpp
  CubitQuadFacetData.hpp
  CurveFacetEvalTool.hpp
  debug.hpp
  FacetDataUtil.hpp
  FacetEntity.hpp
  FacetEvalTool.hpp
  FacetorTool.cpp
  FacetorTool.hpp
  FacetorUtil.hpp
  FacetorUtil.cpp
  GeoNode.hpp
  GeoTet.hpp
  LoopParamTool.hpp
  PointGridSearch.hpp
  PointLoopFacetor.hpp
  TDChordal.hpp
  TDFacetboolData.hpp
  TDFacetBoundaryEdge.hpp
  TDFacetBoundaryPoint.hpp
  TDGeomFacet.hpp
  TDDelaunay.cpp
  TDDelaunay.hpp
  TDInterpNode.hpp
  TetFacetorTool.cpp
  TetFacetorTool.hpp
)

#INCLUDE_DIRECTORIES(${Cholla_SOURCE_DIR})
INCLUDE_DIRECTORIES(${cubit_util_SOURCE_DIR})
INCLUDE_DIRECTORIES(${cubit_util_BINARY_DIR})

ADD_LIBRARY(cholla STATIC ${CHOLLA_SRCS})

install(
  FILES ${CHOLLA_HEADERS}
  DESTINATION include
  COMPONENT Development
)
