
PROJECT(cubit_virtual)

INCLUDE(${cubit_geom_SOURCE_DIR}/../util/UseUtil.cmake)

SET(VIRTUAL_SRCS
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
    SubCurve.cpp
    SubEntitySet.cpp
    SubSurface.cpp
    TDVGFacetOwner.cpp
    TDVGFacetSplit.cpp
    VirtualQueryEngine.cpp
    VirtualImprintTool.cpp
    )

IF(ALPHA_SPLIT_VIRTUAL)
SET(VIRTUAL_SRCS ${VIRTUAL_SRCS} SplitSurfaceVirtual.cpp)
ENDIF(ALPHA_SPLIT_VIRTUAL)

SET(EXTRA_VIRTUAL_SRCS
    AllocMemManagersVirtual.cpp
    )

IF(TEMPLATE_DEFS_INCLUDED)
SET(TEMPLATE_SRCS
    VGArray.cpp
    VGLoopTool.cpp
    )
ENDIF(TEMPLATE_DEFS_INCLUDED)
    
IF(WIN32)
  FOREACH(var ${VIRTUAL_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(VIRTUAL_HEADERS ${VIRTUAL_HEADERS} ${header})
  ENDFOREACH(var ${VIRTUAL_SRCS})
ENDIF(WIN32)

INCLUDE_DIRECTORIES(${cubit_geom_SOURCE_DIR}/facet
                    ${cubit_geom_SOURCE_DIR}/Cholla)

ADD_LIBRARY(cubit_virtual ${VIRTUAL_SRCS} ${EXTRA_VIRTUAL_SRCS} 
                          ${TEMPLATE_SRCS} ${VIRTUAL_HEADERS})

