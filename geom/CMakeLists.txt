

PROJECT(cubit_geom)

OPTION(CGM_FACET "Build CGM with Facet support" ON)
IF(CGM_FACET)
SUBDIRS(Cholla facetbool facet)
ENDIF(CGM_FACET)

OPTION(CGM_VIRTUAL "Build CGM with Virtual support" ON)
IF(CGM_VIRTUAL)
SUBDIRS(virtual)
ENDIF(CGM_VIRTUAL)

# build cgm with acis support, default to off
OPTION(CGM_ACIS "Build CGM with ACIS support" OFF)
IF(CGM_ACIS)
SUBDIRS(ACIS)
ENDIF(CGM_ACIS)


OPTION(CGM_SMLIB "Build CGM with SMLIB support" OFF)
IF(CGM_SMLIB)
SUBDIRS(smlib)
ENDIF(CGM_SMLIB)


# build cgm with granite support, default to off
OPTION(CGM_GRANITE "Build CGM with Granite support" OFF)
IF(CGM_GRANITE)
  SUBDIRS(granite)
ENDIF(CGM_GRANITE)

# build cgm with OpenCASCADE support, default to on
OPTION(CGM_OCC "Build CGM with OpenCASCADE support" ON)
IF(CGM_OCC)
  SUBDIRS(OCC)
ENDIF(CGM_OCC)
OPTION(CGM_CATIA "Build CGM with CATIA support" OFF)
IF(CGM_CATIA)
  SET(CGM_ACIS OFF CACHE BOOL "ACIS not supported with CATIA" FORCE)
  SET(FEATURE_OPTION OFF CACHE BOOL "Feature option not allowed with CATIA" FORCE)
  SET(TOLERANT_TRI OFF CACHE BOOL "Geometry-tolerant meshing option not allowed with CATIA" FORCE)
  SET(ACIS_IGES OFF CACHE BOOL "ACIS_IGES not allowed with CATIA" FORCE)
  SET(ACIS_STEP OFF CACHE BOOL "ACIS_STEP not allowed with CATIA" FORCE)
  SET(CUBIT_GEOM_BUILD_SHARED_LIBS ON CACHE BOOL "GEOM_LIBS need to be shared when CATIA is used" FORCE)
  SET(CUBIT_UTIL_BUILD_SHARED_LIBS ON CACHE BOOL "UTIL_LIBS need to be shared when CATIA is used" FORCE)

ENDIF(CGM_CATIA)

INCLUDE(${cubit_geom_SOURCE_DIR}/../util/UseUtil.cmake)

SET(GEOM_SRCS
    AnalyticGeometryTool.cpp
    AutoMidsurfaceTool.cpp
    BasicTopologyEntity.cpp
    Body.cpp
    BodySM.cpp
    BoundingBoxTool.cpp
    BridgeManager.cpp
    CAActuateSet.cpp
    CADeferredAttrib.cpp
    CAEntityColor.cpp
    CAEntityId.cpp
    CAEntityName.cpp
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

SET(EXTRA_GEOM_SRCS
    AllocMemManagersGeom.cpp
    )


IF(TEMPLATE_DEFS_INCLUDED)
SET(TEMPLATE_SRCS
    BoundaryConstrainTool.cpp
    )
ENDIF(TEMPLATE_DEFS_INCLUDED)

IF(WIN32)
  FOREACH(var ${GEOM_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(GEOM_HEADERS ${GEOM_HEADERS} ${header})
  ENDFOREACH(var ${GEOM_SRCS})
ENDIF(WIN32)

# don't add an include path to other parts of cgm
SET(THIS_INCLUDES 
          ${cubit_util_SOURCE_DIR} 
          ${cubit_util_BINARY_DIR} 
          ${cubit_geom_SOURCE_DIR} 
          ${cubit_geom_BINARY_DIR})

SET_DIRECTORY_PROPERTIES(PROPERTIES INCLUDE_DIRECTORIES "${THIS_INCLUDES}")

OPTION(CUBIT_GEOM_BUILD_SHARED_LIBS "Build cubit_geom as a shared library" OFF)

IF(BUILD_SHARED_LIBS)
  SET(CUBIT_GEOM_BUILD_SHARED_LIBS ON CACHE BOOL "Build cubit_geom as shared library" FORCE)
ENDIF(BUILD_SHARED_LIBS)

IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  # util must be shared if geom is.
  SET(CUBIT_UTIL_BUILD_SHARED_LIBS ON CACHE BOOL "Build cubit util as a shared library" FORCE)
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)


CONFIGURE_FILE(${cubit_geom_SOURCE_DIR}/CubitGeomConfigure.h.in
               ${cubit_geom_BINARY_DIR}/CubitGeomConfigure.h
               @ONLY)

IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_geom SHARED ${GEOM_SRCS} ${EXTRA_GEOM_SRCS} ${GEOM_HEADERS} ${TEMPLATE_SRCS})
ELSE(CUBIT_GEOM_BUILD_SHARED_LIBS)
  ADD_LIBRARY(cubit_geom ${GEOM_SRCS} ${EXTRA_GEOM_SRCS} ${GEOM_HEADERS} ${TEMPLATE_SRCS})
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)

TARGET_LINK_LIBRARIES(cubit_geom cubit_util)

IF(CUBIT_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cubit_geom
                        PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
ENDIF(CUBIT_LIBRARY_PROPERTIES)

SET(CMAKE_INSTALL_BINARY_DIR "bin" CACHE PATH "Install directory for binaries")

IF(CUBIT_GEOM_BUILD_SHARED_LIBS)
  INSTALL(TARGETS cubit_geom
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
          )
ENDIF(CUBIT_GEOM_BUILD_SHARED_LIBS)

IF(BUILD_TESTING)
  #ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

