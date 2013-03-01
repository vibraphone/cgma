# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories
# The values checked are defined in AM_CONDITIONAL statements
# in configure.in
SUBDIRS = . Cholla facetbool facet virtual
if build_ACIS
  SUBDIRS += ACIS
endif
if build_parallel
  SUBDIRS += parallel
endif
if build_OCC
  SUBDIRS += OCC
endif

# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)
AM_CPPFLAGS = -I$(top_srcdir)/util -I$(top_srcdir)/amendment -I.

# The name of the library to build
noinst_LTLIBRARIES =
if BUILD_CGM
  noinst_LTLIBRARIES += libcubit_geom.la
endif

# The directory where headers will be installed
libcubit_geom_la_includedir = $(includedir)

libcubit_geom_la_LIBADD = $(top_builddir)/amendment/libcgma_amendment.la

LDADD = $(top_builddir)/libcgm.la $(CGM_EXT_LIBS) $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS) $(top_builddir)/util/libcubit_util.la

# The non-template sources
if BUILD_CGM
libcubit_geom_la_SOURCES = \
    AllocMemManagersGeom.cpp \
    AnalyticGeometryTool.cpp \
    BasicTopologyEntity.cpp \
    Body.cpp \
    BodySM.cpp \
    BoundingBoxTool.cpp \
    BridgeManager.cpp \
    CAActuateSet.cpp \
    CADeferredAttrib.cpp \
    CAEntityColor.cpp \
    CAEntityId.cpp \
    CAEntityName.cpp \
    CAGroup.cpp \
    CAMergePartner.cpp \
    CAMergeStatus.cpp \
    CASourceFeature.cpp \
    CAUniqueId.cpp \
    CGMApp.cpp \
    CGMEngineDynamicLoader.cpp \
    Chain.cpp \
    CoEdge.cpp \
    CoEdgeSM.cpp \
    CoFace.cpp \
    CollectionEntity.cpp \
    CoVertex.cpp \
    CoVolume.cpp \
    CubitAttrib.cpp \
    CubitAttribManager.cpp \
    CubitAttribUser.cpp \
    CubitPolygon.cpp \
    CubitSimpleAttrib.cpp \
    Curve.cpp \
    CurveOverlapFacet.cpp \
    CurveSM.cpp \
    CylinderEvaluator.cpp \
    DAG.cpp \
    DagDrawingTool.cpp \
    GeomDataObserver.cpp \
    GeometryEntity.cpp \
    GeometryFeatureEngine.cpp \
    GeometryFeatureTool.cpp \
    GeometryHealerEngine.cpp \
    GeometryHealerTool.cpp \
    GeometryModifyEngine.cpp \
    GeometryModifyTool.cpp \
    GeometryQueryEngine.cpp \
    GeometryQueryTool.cpp \
    GeometryUtil.cpp \
    GeomMeasureTool.cpp \
    GroupingEntity.cpp \
    GSaveOpen.cpp \
    Loop.cpp \
    LoopSM.cpp \
    Lump.cpp \
    LumpSM.cpp \
    MedialTool2D.cpp \
    MedialTool3D.cpp \
    MergeToolAssistant.cpp \
    MergeTool.cpp \
    MidPlaneTool.cpp \
    ModelEntity.cpp \
    ModelQueryEngine.cpp \
    OffsetSplitTool.cpp \
    OldUnmergeCode.cpp \
    Point.cpp \
    PointSM.cpp \
    RefCollection.cpp \
    RefEdge.cpp \
    RefEntity.cpp \
    RefEntityFactory.cpp \
    RefEntityName.cpp \
    RefFace.cpp \
    RefGroup.cpp \
    RefVertex.cpp \
    RefVolume.cpp \
    SenseEntity.cpp \
    Shell.cpp \
    ShellSM.cpp \
    SphereEvaluator.cpp \
    SplitSurfaceTool.cpp \
    Surface.cpp \
    SurfaceOverlapFacet.cpp \
    SurfaceOverlapTool.cpp \
    SurfParamTool.cpp \
    SurfaceSM.cpp \
    TDCAGE.cpp \
    TDSourceFeature.cpp \
    TDSplitSurface.cpp \
    TDSurfaceOverlap.cpp \
    TDUniqueId.cpp \
    TopologyBridge.cpp \
    TopologyEntity.cpp 
endif

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_geom_la_include_HEADERS = \
  AnalyticGeometryTool.hpp \
  BasicTopologyEntity.hpp \
  Body.hpp \
  BodySM.hpp \
  BoundingBoxTool.hpp \
  BridgeManager.hpp \
  CAActuateSet.hpp \
  CADeferredAttrib.hpp \
  CADefines.hpp \
  CAEntityColor.hpp \
  CAEntityId.hpp \
  CAEntityName.hpp \
  CAGroup.hpp \
  CAMergePartner.hpp \
  CAMergeStatus.hpp \
  CASourceFeature.hpp \
  CAUniqueId.hpp \
  CGMApp.hpp \
  CGMEngineDynamicLoader.hpp \
  Chain.hpp \
  CoEdge.hpp \
  CoEdgeSM.hpp \
  CoFace.hpp \
  CoVertex.hpp \
  CoVolume.hpp \
  CollectionEntity.hpp \
  CompositeCombineEvent.hpp \
  CubitAttrib.hpp \
  CubitAttribManager.hpp \
  CubitAttribUser.hpp \
  CubitEvaluator.hpp \
  CubitGeomConfigure.h \
  CubitPolygon.hpp \
  CubitSimpleAttrib.hpp \
  Curve.hpp \
  CurveOverlapFacet.hpp \
  CurveSM.hpp \
  CylinderEvaluator.hpp \
  DAG.hpp \
  DagDrawingTool.hpp \
  DagType.hpp \
  GSaveOpen.hpp \
  GeomDataObserver.hpp \
  GeomMeasureTool.hpp \
  GeomPoint.hpp \
  GeomSeg.hpp \
  GeometryEntity.hpp \
  GeometryFeatureEngine.hpp \
  GeometryFeatureTool.hpp \
  GeometryHealerEngine.hpp \
  GeometryHealerTool.hpp \
  GeometryModifyEngine.hpp \
  GeometryModifyTool.hpp \
  GeometryQueryEngine.hpp \
  GeometryQueryTool.hpp \
  GeometryUtil.hpp \
  GroupingEntity.hpp \
  IntermediateGeomEngine.hpp \
  Loop.hpp \
  LoopSM.hpp \
  Lump.hpp \
  LumpSM.hpp \
  MedialTool2D.hpp \
  MedialTool3D.hpp \
  MergeTool.hpp \
  MergeToolAssistant.hpp \
  MidPlaneTool.hpp \
  ModelEntity.hpp \
  ModelQueryEngine.hpp \
  OffsetSplitTool.hpp \
  OldUnmergeCode.hpp \
  Point.hpp \
  PointSM.hpp \
  RefCollection.hpp \
  RefEdge.hpp \
  RefEntity.hpp \
  RefEntityFactory.hpp \
  RefEntityName.hpp \
  RefEntityNameMap.hpp \
  RefFace.hpp \
  RefGroup.hpp \
  RefVertex.hpp \
  RefVolume.hpp \
  SenseEntity.hpp \
  Shell.hpp \
  ShellSM.hpp \
  SphereEvaluator.hpp \
  SplitSurfaceTool.hpp \
  SurfParamTool.hpp \
  Surface.hpp \
  SurfaceOverlapFacet.hpp \
  SurfaceOverlapTool.hpp \
  SurfaceSM.hpp \
  TBOwner.hpp \
  TBOwnerSet.hpp \
  TDCAGE.hpp \
  TDCompare.hpp \
  TDSourceFeature.hpp \
  TDSplitSurface.hpp \
  TDSurfaceOverlap.hpp \
  TDUniqueId.hpp \
  TopologyBridge.hpp \
  TopologyEntity.hpp \
  UnMergeEvent.hpp 


