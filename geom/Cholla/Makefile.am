# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)
INCLUDES = -I$(top_srcdir)/util -I$(top_srcdir)/geom

# The name of the library to build
lib_LTLIBRARIES = libcholla.la

# The directory where headers will be installed
libcholla_la_includedir = $(includedir)

# The non-template sources
libcholla_la_SOURCES = \
    AllocMemManagersCholla.cpp \
    Cholla.cpp \
    ChollaCurve.cpp \
    ChollaEngine.cpp \
    ChollaEntity.cpp \
    ChollaPoint.cpp \
    ChollaSkinTool.cpp \
    ChollaSurface.cpp \
    ChordalAxis.cpp \
    CubitFacet.cpp \
    CubitFacetData.cpp \
    CubitFacetEdge.cpp \
    CubitFacetEdgeData.cpp \
    CubitPoint.cpp \
    CubitPointData.cpp \
    CubitQuadFacet.cpp \
    CubitQuadFacetData.cpp \
    CurveFacetEvalTool.cpp \
    debug.cpp \
    FacetDataUtil.cpp \
    FacetEntity.cpp \
    FacetEvalTool.cpp \
    GeoNode.cpp \
    GeoTet.cpp \
    LoopParamTool.cpp \
    PointGridSearch.cpp \
    PointLoopFacetor.cpp \
    TDChordal.cpp \
    TDFacetboolData.cpp \
    TDFacetBoundaryEdge.cpp \
    TDFacetBoundaryPoint.cpp \
    TDGeomFacet.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcholla_la_include_HEADERS = \
    BoundaryConstrainTool.hpp \
    Cholla.h \
    ChollaCurve.hpp \
    ChollaEngine.hpp \
    ChollaEntity.hpp \
    ChollaPoint.hpp \
    ChollaSkinTool.hpp \
    ChollaSurface.hpp \
    ChordalAxis.hpp \
    CubitFacet.hpp \
    CubitFacetData.hpp \
    CubitFacetEdge.hpp \
    CubitFacetEdgeData.hpp \
    CubitPoint.hpp \
    CubitPointData.hpp \
    CubitQuadFacet.hpp \
    CubitQuadFacetData.hpp \
    CurveFacetEvalTool.hpp \
    FacetDataUtil.hpp \
    FacetEntity.hpp \
    FacetEvalTool.hpp \
    FacetorTool.hpp \
    GeoNode.hpp \
    GeoTet.hpp \
    LoopParamTool.hpp \
    PointGridSearch.hpp \
    PointLoopFacetor.hpp \
    TDChordal.hpp \
    TDDelaunay.hpp \
    TDFacetBoundaryEdge.hpp \
    TDFacetBoundaryPoint.hpp \
    TDFacetboolData.hpp \
    TDGeomFacet.hpp \
    TDInterpNode.hpp \
    TetFacetorTool.hpp \
    debug.hpp 

# If template defs are included, then the template definitions
# need to be installed with the headers.  Otherwise they need
# to be compiled.
if INCLUDE_TEMPLATE_DEFS
  libcholla_la_include_HEADERS += \
    BoundaryConstrainTool.cpp \
    FacetorTool.cpp \
    TetFacetorTool.cpp \
    TDDelaunay.cpp 
else
  libcholla_la_SOURCES += \
    BoundaryConstrainTool.cpp \
    FacetorTool.cpp \
    TetFacetorTool.cpp \
    TDDelaunay.cpp 
endif

