project(cubit_geom)

option(CGM_FACET "Build CGM with Facet support" ON)
option(CGM_VIRTUAL "Build CGM with Virtual support" ON)
# build cgm with acis support, default to off
option(CGM_ACIS "Build CGM with ACIS support" OFF)

#turn CGM_KCM on by default on windows systems
if(CMAKE_SYSTEM MATCHES "Windows")
  option(CGM_KCM "Build and use the KCM modeler capabilities" ON)
else()
  option(CGM_KCM "Build and use the KCM modeler capabilities" OFF)
endif()

option(CGM_SMLIB "Build CGM with SMLIB support" OFF)

option(CGM_OCC "Build CGM with Open Cascade support" OFF)

# build cgm with granite support, default to off
if(CAT_BUILD)
  option(CGM_GRANITE "Build CGM with Granite support" ON)
else(CAT_BUILD)
  option(CGM_GRANITE "Build CGM with Granite support" OFF)
endif(CAT_BUILD)

option(CGM_CATIA "Build CGM with CATIA support" OFF)
if(CGM_CATIA)
  set(BUILD_SHARED_LIBS ON CACHE BOOL "GEOM_LIBS need to be shared when CATIA is used" FORCE)
endif(CGM_CATIA)

# all parts of CGM need util
include_directories(
  ${cubit_util_SOURCE_DIR}
  ${cubit_util_BINARY_DIR}
)

if(CGM_FACET)
  add_subdirectory(Cholla)
endif()

# all remaining parts need geom
include_directories(
  ${cubit_geom_SOURCE_DIR}
  ${cubit_geom_BINARY_DIR}
)

set(GEOM_SRCS
    AnalyticGeometryTool.cpp
    AutoMidsurfaceTool.cpp
    BasicTopologyEntity.cpp
    Body.cpp
    BodySM.cpp
    BoundingBoxTool.cpp
    BridgeManager.cpp
    CAActuateSet.cpp
    CAEntitySense.cpp
    CADeferredAttrib.cpp
    CAEntityColor.cpp
    CAEntityId.cpp
    CAEntityName.cpp
    CAEntityTol.cpp
    CAGroup.cpp
    CAMergePartner.cpp
    CAMergeStatus.cpp
    CASourceFeature.cpp
    CAUniqueId.cpp
    CGMApp.cpp
    CGMEngineDynamicLoader.cpp
    CGMHistory.cpp
    Chain.cpp
    CoEdge.cpp
    CoEdgeSM.cpp
    CoFace.cpp
    CollectionEntity.cpp
    CoVertex.cpp
    CoVolume.cpp
    CubitAttrib.cpp
    CubitAttribManager.cpp
    CubitAttribUser.cpp
    CubitCompat.cpp
    CubitPolygon.cpp
    CubitSimpleAttrib.cpp
    Curve.cpp
    CurveOverlapFacet.cpp
    CurveSM.cpp
    CylinderEvaluator.cpp
    DAG.cpp
    DagDrawingTool.cpp
    GeomDataObserver.cpp
    GeometryEntity.cpp
    GeometryFeatureTool.cpp
    GeometryHealerTool.cpp
    GeometryModifyEngine.cpp
    GeometryModifyTool.cpp
    GeometryQueryEngine.cpp
    GeometryQueryTool.cpp
    GeometryUtil.cpp
    GeomMeasureTool.cpp
    GfxPreview.cpp
    GroupingEntity.cpp
    GSaveOpen.cpp
    LocalToleranceTool.cpp
    Loop.cpp
    LoopSM.cpp
    Lump.cpp
    LumpSM.cpp
    MedialTool2D.cpp
    MedialTool3D.cpp
    MergeToolAssistant.cpp
    MergeTool.cpp
    ModelEntity.cpp
    ModelQueryEngine.cpp
    OffsetSplitTool.cpp
    OldUnmergeCode.cpp
    Point.cpp
    PointSM.cpp
    RefCollection.cpp
    RefEdge.cpp
    RefEntity.cpp
    RefEntityFactory.cpp
    RefEntityName.cpp
    RefFace.cpp
    RefGroup.cpp
    RefVertex.cpp
    RefVolume.cpp
    SenseEntity.cpp
    Shell.cpp
    ShellSM.cpp
    SphereEvaluator.cpp
    SplitSurfaceTool.cpp
    Surface.cpp
    SurfaceOverlapFacet.cpp
    SurfaceOverlapTool.cpp
    SurfParamTool.cpp
    SurfaceSM.cpp
    TDCAGE.cpp
    TDSourceFeature.cpp
    TDSplitSurface.cpp
    TDSurfaceOverlap.cpp
    TDUniqueId.cpp
    TopologyBridge.cpp
    TopologyEntity.cpp
    )

if(CAT_BUILD)
  set(GEOM_SRCS ${GEOM_SRCS} CAProWeld.cpp TDProWeld.cpp)
endif(CAT_BUILD)

set(EXTRA_GEOM_SRCS
    AllocMemManagersGeom.cpp
    )

# Create a list of header-only files:
set(GEOM_HEADERS
  CADefines.hpp
  CompositeCombineEvent.hpp
  CubitCompat.hpp
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h" # Generated
  CubitEvaluator.hpp
  CubitGeomConfigure.h
  DagType.hpp
  GeomPoint.hpp
  GeomSeg.hpp
  GeometryFeatureEngine.hpp
  GeometryHealerEngine.hpp
  IntermediateGeomEngine.hpp
  MidPlaneTool.hpp
  RefEntityNameMap.hpp
  TBOwner.hpp
  TBOwnerSet.hpp
  TDCompare.hpp
  UnMergeEvent.hpp
)
# Now algorithmically generate header names
# from implementation filenames:
foreach(var ${GEOM_SRCS})
  string(REGEX REPLACE ".cpp" ".hpp" header ${var})
  set(GEOM_HEADERS ${GEOM_HEADERS} ${header})
endforeach(var ${GEOM_SRCS})

add_library(cubit_geom ${GEOM_SRCS} ${EXTRA_GEOM_SRCS})
target_link_libraries(cubit_geom cubit_util)
set(CGMA_LIBS "${CGMA_LIBS}" cubit_geom)

if(CUBIT_LIBRARY_PROPERTIES)
  set_target_properties(cubit_geom
    PROPERTIES ${CUBIT_LIBRARY_PROPERTIES}
  )
endif(CUBIT_LIBRARY_PROPERTIES)

set(CMAKE_INSTALL_BINARY_DIR "lib" CACHE PATH "Install directory for binaries")
set(CMAKE_INSTALL_INCLUDE_DIR "include" CACHE PATH "Install directory for headers")

install(
  TARGETS cubit_geom  ${CUBIT_GEOM_EXPORT_GROUP}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)

if(CGM_FACET)
  add_subdirectory(facetbool)
  add_subdirectory(facet)
endif(CGM_FACET)

if(CGM_VIRTUAL)
  add_subdirectory(virtual)
endif(CGM_VIRTUAL)

if(CGM_ACIS)
  add_subdirectory(ACIS)
endif(CGM_ACIS)

if(CGM_SMLIB)
  add_subdirectory(smlib)
endif(CGM_SMLIB)

if(CGM_OCC)
  add_subdirectory(OCC)
endif(CGM_OCC)

if(CGM_GRANITE)
  add_subdirectory(granite)
endif(CGM_GRANITE)

if(CGM_CATIA)
  add_subdirectory(catia)
endif(CGM_CATIA)

if(CGM_KCM)
  add_subdirectory(KCM)
endif(CGM_KCM)

configure_file(CGMConfigure.h.in "${CMAKE_CURRENT_BINARY_DIR}/CGMConfigure.h")

set(CUBIT_GEOM_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  ${cubit_geom_SOURCE_DIR}/CubitGeomConfigure.h.in
  ${cubit_geom_BINARY_DIR}/CubitGeomConfigure.h
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CubitCompat.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h"
  @ONLY
)

# Pass to parent environment
set(CGMA_LIBS "${CGMA_LIBS}" PARENT_SCOPE)
set(CGMA_DEPLIBS ${CGMA_DEPLIBS} ${CGMA_OCC_DEPLIBS} ${CGMA_ACIS_DEPLIBS} PARENT_SCOPE)
set(CGMA_ACIS_DEFINES "${CGMA_ACIS_DEFINES}" PARENT_SCOPE)
set(CGMA_OCC_DEFINES "${CGMA_OCC_DEFINES}" PARENT_SCOPE)
set(OCE_INCLUDE_DIRS ${OCE_INCLUDE_DIRS} PARENT_SCOPE)

add_definitions(${CGMA_ACIS_DEFINES} ${CGMA_OCC_DEFINES})

install(
  FILES ${GEOM_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDE_DIR}
  COMPONENT Development)

#x3d directory
find_path(X3D_DIR "" DOC "X3D directory.")


if(BUILD_TESTING)
  #add_subdirectory(testing)
endif(BUILD_TESTING)
