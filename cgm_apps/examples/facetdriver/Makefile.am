# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_ACIS_DEF) $(CUBIT_CGM)

TESTS = facets
EXTRA_DIST = README

check_PROGRAMS = facets

facets_SOURCES = facets.cpp

AM_CPPFLAGS += -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/Cholla \
               -I$(top_srcdir)/geom/facet

LDADD = $(top_builddir)/libcgm.la $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
