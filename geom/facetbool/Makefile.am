# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)
AM_CPPFLAGS = -I$(top_srcdir)/util \
           -I$(top_srcdir)/amendment \
           -I$(top_srcdir)/geom 

# The name of the library to build
lib_LTLIBRARIES = 
if BUILD_CGM
  lib_LTLIBRARIES += libcubit_facetbool.la
endif

# The directory where headers will be installed
libcubit_facetbool_la_includedir = $(includedir)

# The non-template sources
libcubit_facetbool_la_SOURCES = \
    FBClassify.cpp \
    FBDataUtil.cpp \
    FBImprint.cpp \
    FBIntersect.cpp \
    FBPolyhedron.cpp \
    FBRetriangulate.cpp \
    FBTiler.cpp \
    IntegerHash.cpp \
    KdTree.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_facetbool_la_include_HEADERS = \
    FBClassify.hpp \
    FBDataUtil.hpp \
    FBDefines.hpp \
    FBImprint.hpp \
    FBIntersect.hpp \
    FBPolyhedron.hpp \
    FBRetriangulate.hpp \
    FBStructs.hpp \
    FBTiler.hpp \
    IntegerHash.hpp \
    KdTree.hpp 

