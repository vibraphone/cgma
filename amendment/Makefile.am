# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(top_srcdir)/geom -I$(top_srcdir)/util -I.

noinst_LTLIBRARIES = libcgma_amendment.la
libcgma_amendment_la_includedir = $(includedir)
libcgma_amendment_la_SOURCES = AutoMidsurfaceTool.cpp \
    CAEntityTol.cpp \
    CAEntitySense.cpp \
    CGMHistory.cpp \
    CubitCompat.cpp \
    GfxPreview.cpp \
    LocalToleranceTool.cpp \
    CubitFileSimModel.cpp \
    CubitUndo.cpp

libcgma_amendment_la_include_HEADERS = AutoMidsurfaceTool.hpp \
    CAEntityTol.hpp \
    CAEntitySense.cpp \
    CGMHistory.hpp \
    CubitCompat.hpp \
    CubitCompat.h \
    GfxPreview.hpp \
    LocalToleranceTool.hpp \
    CubitFileSimModel.hpp \
    CubitUndo.hpp

if BUILD_CGM
  LDADD = $(top_builddir)/libcgm.la $(CGM_EXT_LIBS) $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS) 
else
  LDADD = -lcubit_util -lcubit_geom -lcubiti19 $(CGM_EXT_LIBS) $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
endif


