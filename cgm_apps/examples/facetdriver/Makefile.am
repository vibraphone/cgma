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
               -I$(top_srcdir)/geom/facetbool

LDADD = $(top_builddir)/libcgm.la $(CGM_EXT_LIBS) 
AM_LDFLAGS += $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)  -lcubit_geom -lcubit_util -lcubit_facetbool -lcubit_facet -lCholla -lcubit_virtual -ldl
