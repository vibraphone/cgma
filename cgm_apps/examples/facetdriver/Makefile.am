#DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_ACIS_DEF) $(CUBIT_CGM)
TESTS = facets
EXTRA_DIST = README

check_PROGRAMS = facets

facets_SOURCES = facets.cpp

AM_CPPFLAGS += $(HAVE_ACIS_DEF) $(CUBIT_CGM) \
               -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/Cholla \
               -I$(top_srcdir)/geom/facet \
               -I$(top_srcdir)/geom/facetbool \
               -I$(top_srcdir)/amendment 

if BUILD_CGM
  LDADD = $(top_builddir)/libcgm.la $(CGM_EXT_LIBS) $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS) $(top_builddir)/geom/Cholla/libCholla.la $(top_builddir)/geom/facet/libcubit_facet.la $(top_builddir)/geom/facetbool/libcubit_facetbool.la $(top_builddir)/geom/libcubit_geom.la $(top_builddir)/util/libcubit_util.la 
else
  LDADD = -lcubit_geom -lcubit_util -lcubiti19 $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
endif
