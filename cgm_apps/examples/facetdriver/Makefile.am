TESTS = facets_90
EXTRA_DIST = README $(TESTS)

check_PROGRAMS = facets

facets_SOURCES = facets.cpp

AM_CPPFLAGS += -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/Cholla \
               -I$(top_srcdir)/geom/facet

LDADD = ../../../geom/libcubit_geom.la \
        ../../../geom/facet/libcubit_facet.la 
