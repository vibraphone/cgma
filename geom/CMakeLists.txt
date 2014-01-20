PROJECT(cubit_geom)

OPTION(CGM_FACET "Build CGM with Facet support" ON)
OPTION(CGM_VIRTUAL "Build CGM with Virtual support" ON)
# build cgm with acis support, default to off
OPTION(CGM_ACIS "Build CGM with ACIS support" OFF)

#turn CGM_KCM on by default on windows systems
if(CMAKE_SYSTEM MATCHES "Windows")
  OPTION(CGM_KCM "Build and use the KCM modeler capabilities" ON)
else()
  OPTION(CGM_KCM "Build and use the KCM modeler capabilities" OFF)
endif()

OPTION(CGM_SMLIB "Build CGM with SMLIB support" OFF)

OPTION(CGM_OCC "Build CGM with Open Cascade support" OFF)

# build cgm with granite support, default to off
IF(CAT_BUILD)
  OPTION(CGM_GRANITE "Build CGM with Granite support" ON)
ELSE(CAT_BUILD)
  OPTION(CGM_GRANITE "Build CGM with Granite support" OFF)
ENDIF(CAT_BUILD)

OPTION(CGM_CATIA "Build CGM with CATIA support" OFF)
IF(CGM_CATIA)
  SET(BUILD_SHARED_LIBS ON CACHE BOOL "GEOM_LIBS need to be shared when CATIA is used" FORCE)
ENDIF(CGM_CATIA)



# all parts of CGM need util
include_directories(
          ${cubit_util_SOURCE_DIR} 
          ${cubit_util_BINARY_DIR} 
          )

IF(CGM_FACET)
  ADD_SUBDIRECTORY(Cholla)
ENDIF()

# all remaining parts need geom
include_directories(
          ${cubit_geom_SOURCE_DIR} 
          ${cubit_geom_BINARY_DIR} 
          )

SET(GEOM_SRCS
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
    PeriodicParamTool.cpp
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

IF(CAT_BUILD)
  SET(GEOM_SRCS ${GEOM_SRCS} CAProWeld.cpp TDProWeld.cpp)
ENDIF(CAT_BUILD)

SET(EXTRA_GEOM_SRCS
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
FOREACH(var ${GEOM_SRCS})
  STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
  SET(GEOM_HEADERS ${GEOM_HEADERS} ${header})
ENDFOREACH(var ${GEOM_SRCS})

ADD_LIBRARY(cubit_geom ${GEOM_SRCS} ${EXTRA_GEOM_SRCS})
TARGET_LINK_LIBRARIES(cubit_geom cubit_util)

IF(CUBIT_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cubit_geom
                        PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
ENDIF(CUBIT_LIBRARY_PROPERTIES)

SET(CMAKE_INSTALL_BINARY_DIR "bin" CACHE PATH "Install directory for binaries")
SET(CMAKE_INSTALL_INCLUDE_DIR "include" CACHE PATH "Install directory for headers")

IF(BUILD_SHARED_LIBS)
  INSTALL(TARGETS cubit_geom  ${CUBIT_GEOM_EXPORT_GROUP}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
          )
ENDIF(BUILD_SHARED_LIBS)


IF(CGM_FACET)
  ADD_SUBDIRECTORY(facetbool)
  ADD_SUBDIRECTORY(facet)
ENDIF(CGM_FACET)

IF(CGM_VIRTUAL)
  ADD_SUBDIRECTORY(virtual)
ENDIF(CGM_VIRTUAL)

IF(CGM_ACIS)
  ADD_SUBDIRECTORY(ACIS)
ENDIF(CGM_ACIS)

IF(CGM_SMLIB)
  ADD_SUBDIRECTORY(smlib)
ENDIF(CGM_SMLIB)

IF(CGM_OCC)
  ADD_SUBDIRECTORY(OCC)
ENDIF(CGM_OCC)

IF(CGM_GRANITE)
  ADD_SUBDIRECTORY(granite)
ENDIF(CGM_GRANITE) 

IF(CGM_CATIA)
  ADD_SUBDIRECTORY(catia)
ENDIF(CGM_CATIA)

IF(CGM_KCM)
  ADD_SUBDIRECTORY(KCM)
ENDIF(CGM_KCM)

CONFIGURE_FILE(CGMConfigure.h.in "${CMAKE_CURRENT_BINARY_DIR}/CGMConfigure.h")

SET(CUBIT_GEOM_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
CONFIGURE_FILE(${cubit_geom_SOURCE_DIR}/CubitGeomConfigure.h.in
               ${cubit_geom_BINARY_DIR}/CubitGeomConfigure.h
               @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CubitCompat.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h"
  @ONLY
)

install(
  FILES ${GEOM_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDE_DIR}
  COMPONENT Development)

#x3d directory
FIND_PATH(X3D_DIR "" DOC "X3D directory.")


IF(BUILD_TESTING)
  #ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
