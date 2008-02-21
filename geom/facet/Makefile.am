# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)
INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           -I$(top_srcdir)/geom/Cholla \
           -I$(top_srcdir)/geom/facetbool 

# The name of the library to build
lib_LTLIBRARIES = 
if BUILD_CGM
  lib_LTLIBRARIES += libcubit_facet.la
endif

libcubit_facet_la_LIBADD = ../libcubit_geom.la \
                           ../facetbool/libcubit_facetbool.la \
                           ../Cholla/libcholla.la 

# The directory where headers will be installed
libcubit_facet_la_includedir = $(includedir)

# The non-template sources
libcubit_facet_la_SOURCES = \
    FacetAttrib.cpp \
    FacetAttribSet.cpp \
    FacetBody.cpp \
    FacetboolInterface.cpp \
    FacetCoEdge.cpp \
    FacetCurve.cpp \
    FacetLoop.cpp \
    FacetLump.cpp \
    FacetModifyEngine.cpp \
    FacetParamTool.cpp \
    FacetPoint.cpp \
    FacetQueryEngine.cpp \
    FacetShell.cpp \
    FacetSurface.cpp \
    GridSearchTree.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_facet_la_include_HEADERS = \
    FacetAttrib.hpp \
    FacetAttribSet.hpp \
    FacetBody.hpp \
    FacetCoEdge.hpp \
    FacetCurve.hpp \
    FacetGeometryCreator.hpp \
    FacetLoop.hpp \
    FacetLump.hpp \
    FacetModifyEngine.hpp \
    FacetParamTool.hpp \
    FacetPoint.hpp \
    FacetQueryEngine.hpp \
    FacetShell.hpp \
    FacetSurface.hpp \
    FacetboolInterface.hpp \
    GridSearchTree.hpp \
    GridSearchTreeNode.hpp 

