EXTRA_DIST = README

check_PROGRAMS = facets

facets_SOURCES = facets.cpp

TESTS = facets_90

CPPFLAGS += -I$(top_srcdir)/util \
            -I$(top_srcdir)/geom \
            -I$(top_srcdir)/geom/Cholla \
            -I$(top_srcdir)/geom/facet

LDADD = ../../../geom/libcubit_geom.la \
        ../../../geom/facet/libcubit_facet.la 
