project(cubit_virtual)

set(VIRTUAL_SRCS
  CACompositeVG.cpp
  CAPartitionVG.cpp
  CAVirtualVG.cpp
  CollapseAngleTool.cpp
  CollapseCurveTool.cpp
  CompositeAttrib.cpp
  CompositeBody.cpp
  CompositeCoEdge.cpp
  CompositeCoSurf.cpp
  CompositeCurve.cpp
  CompositeEngine.cpp
  CompositeGeom.cpp
  CompositeLoop.cpp
  CompositeLump.cpp
  CompositePoint.cpp
  CompositeShell.cpp
  CompositeSurface.cpp
  CompositeTool.cpp
  CompSurfFacets.cpp
  Faceter.cpp
  FaceterFacetData.cpp
  FaceterPointData.cpp
  FacetProjectTool.cpp
  HiddenEntitySet.cpp
  ImprintBoundaryTool.cpp
  ImprintLineSegment.cpp
  ImprintPointData.cpp
  PartitionBody.cpp
  PartitionCoEdge.cpp
  PartitionCoSurf.cpp
  PartitionCurve.cpp
  PartitionEngine.cpp
  PartitionEntity.cpp
  PartitionLoop.cpp
  PartitionLump.cpp
  PartitionLumpImprint.cpp
  PartitionPoint.cpp
  PartitionShell.cpp
  PartitionSurface.cpp
  PartitionTool.cpp
  PartPTCurve.cpp
  PartSurfFacetTool.cpp
  PST_Data.cpp
  RemoveBlends.cpp
  SegmentedCurve.cpp
  SimplifyTool.cpp
  SplitCompositeSurfaceTool.cpp
  SubCurve.cpp
  SubEntitySet.cpp
  SubSurface.cpp
  TDVGFacetOwner.cpp
  TDVGFacetSplit.cpp
  VirtualQueryEngine.cpp
  VirtualImprintTool.cpp
)

if(ALPHA_SPLIT_VIRTUAL)
  set(VIRTUAL_SRCS ${VIRTUAL_SRCS} SplitSurfaceVirtual.cpp)
endif()

set(EXTRA_VIRTUAL_SRCS
  AllocMemManagersVirtual.cpp
)

set(VIRTUAL_HEADERS
  CACompositeVG.hpp
  CAPartitionVG.hpp
  CAVirtualVG.hpp
  CollapseAngleTool.hpp
  CompSurfFacets.hpp
  CompositeAttrib.hpp
  CompositeBody.hpp
  CompositeCoEdge.hpp
  CompositeCoSurf.hpp
  CompositeCurve.hpp
  CompositeEngine.hpp
  CompositeGeom.hpp
  CompositeLoop.hpp
  CompositeLump.hpp
  CompositePoint.hpp
  CompositeShell.hpp
  CompositeSurface.hpp
  CompositeTool.hpp
  FacetProjectTool.hpp
  Faceter.hpp
  FaceterFacetData.hpp
  FaceterPointData.hpp
  HiddenEntitySet.hpp
  ImprintBoundaryTool.hpp
  ImprintLineSegment.hpp
  ImprintMatchData.hpp
  ImprintPointData.hpp
  PST_Data.hpp
  PartPTCurve.hpp
  PartSurfFacetTool.hpp
  PartitionBody.hpp
  PartitionCoEdge.hpp
  PartitionCoSurf.hpp
  PartitionCurve.hpp
  PartitionEngine.hpp
  PartitionEntity.hpp
  PartitionLoop.hpp
  PartitionLump.hpp
  PartitionLumpImprint.hpp
  PartitionPoint.hpp
  PartitionShell.hpp
  PartitionSurface.hpp
  PartitionTool.hpp
  SegmentedCurve.hpp
  SplitSurfaceVirtual.hpp
  SubCurve.hpp
  SubEntitySet.hpp
  SubSurface.hpp
  TDVGFacetOwner.hpp
  TDVGFacetSplit.hpp
  VGArray.cpp
  VGArray.hpp
  VGDefines.h
  VGLoopTool.cpp
  VGLoopTool.hpp
  VirtualImprintTool.hpp
  VirtualQueryEngine.hpp
)

include_directories(
  ${cubit_geom_SOURCE_DIR}/facet
  ${cubit_geom_SOURCE_DIR}/Cholla)

add_library(cubit_virtual STATIC ${VIRTUAL_SRCS} ${EXTRA_VIRTUAL_SRCS})

install(
  TARGETS cubit_virtual
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)
install(
  FILES ${VIRTUAL_HEADERS}
  DESTINATION include
  COMPONENT Development
)
