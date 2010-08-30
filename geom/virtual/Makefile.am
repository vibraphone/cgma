AM_CPPFLAGS += -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/facet \
               -I$(top_srcdir)/geom/Cholla

# The name of the library to build
noinst_LTLIBRARIES = 
if BUILD_CGM
  noinst_LTLIBRARIES += libcubit_virtual.la
endif

# The directory where headers will be installed
libcubit_virtual_la_includedir = $(includedir)

# The non-template sources
libcubit_virtual_la_SOURCES = \
    AllocMemManagersVirtual.cpp \
    CACompositeVG.cpp \
    CAPartitionVG.cpp \
    CAVirtualVG.cpp \
    CollapseAngleTool.cpp \
    CompositeAttrib.cpp \
    CompositeBody.cpp \
    CompositeCoEdge.cpp \
    CompositeCoSurf.cpp \
    CompositeCurve.cpp \
    CompositeEngine.cpp \
    CompositeGeom.cpp \
    CompositeLoop.cpp \
    CompositeLump.cpp \
    CompositePoint.cpp \
    CompositeShell.cpp \
    CompositeSurface.cpp \
    CompositeTool.cpp \
    CompSurfFacets.cpp \
    Faceter.cpp \
    FaceterFacetData.cpp \
    FaceterPointData.cpp \
    FacetProjectTool.cpp \
    HiddenEntitySet.cpp \
    ImprintBoundaryTool.cpp \
    ImprintLineSegment.cpp \
    ImprintPointData.cpp \
    PartitionBody.cpp \
    PartitionCoEdge.cpp \
    PartitionCoSurf.cpp \
    PartitionCurve.cpp \
    PartitionEngine.cpp \
    PartitionEntity.cpp \
    PartitionLoop.cpp \
    PartitionLump.cpp \
    PartitionLumpImprint.cpp \
    PartitionPoint.cpp \
    PartitionShell.cpp \
    PartitionSurface.cpp \
    PartitionTool.cpp \
    PartPTCurve.cpp \
    PartSurfFacetTool.cpp \
    PST_Data.cpp \
    SegmentedCurve.cpp \
    SubCurve.cpp \
    SubEntitySet.cpp \
    SubSurface.cpp \
    TDVGFacetOwner.cpp \
    TDVGFacetSplit.cpp \
    VirtualQueryEngine.cpp \
    VirtualImprintTool.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_virtual_la_include_HEADERS = \
    CACompositeVG.hpp \
    CAPartitionVG.hpp \
    CAVirtualVG.hpp \
    CollapseAngleTool.hpp \
    CompSurfFacets.hpp \
    CompositeAttrib.hpp \
    CompositeBody.hpp \
    CompositeCoEdge.hpp \
    CompositeCoSurf.hpp \
    CompositeCurve.hpp \
    CompositeEngine.hpp \
    CompositeGeom.hpp \
    CompositeLoop.hpp \
    CompositeLump.hpp \
    CompositePoint.hpp \
    CompositeShell.hpp \
    CompositeSurface.hpp \
    CompositeTool.hpp \
    FacetProjectTool.hpp \
    Faceter.hpp \
    FaceterFacetData.hpp \
    FaceterPointData.hpp \
    HiddenEntitySet.hpp \
    ImprintBoundaryTool.hpp \
    ImprintLineSegment.hpp \
    ImprintMatchData.hpp \
    ImprintPointData.hpp \
    PST_Data.hpp \
    PartPTCurve.hpp \
    PartSurfFacetTool.hpp \
    PartitionBody.hpp \
    PartitionCoEdge.hpp \
    PartitionCoSurf.hpp \
    PartitionCurve.hpp \
    PartitionEngine.hpp \
    PartitionEntity.hpp \
    PartitionLoop.hpp \
    PartitionLump.hpp \
    PartitionLumpImprint.hpp \
    PartitionPoint.hpp \
    PartitionShell.hpp \
    PartitionSurface.hpp \
    PartitionTool.hpp \
    SegmentedCurve.hpp \
    SplitSurfaceVirtual.hpp \
    SubCurve.hpp \
    SubEntitySet.hpp \
    SubSurface.hpp \
    TDVGFacetOwner.hpp \
    TDVGFacetSplit.hpp \
    VGArray.hpp \
    VGDefines.h \
    VGLoopTool.hpp \
    VirtualImprintTool.hpp \
    VirtualQueryEngine.hpp 

# If template defs are included, then the template definitions
# need to be installed with the headers.  Otherwise they need
# to be compiled.
if INCLUDE_TEMPLATE_DEFS
  libcubit_virtual_la_include_HEADERS += \
    VGArray.cpp \
    VGLoopTool.cpp 
else
  libcubit_virtual_la_SOURCES += \
    VGArray.cpp \
    VGLoopTool.cpp 
endif

