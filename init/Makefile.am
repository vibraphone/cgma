# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign


# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_ACIS_DEF) $(CUBIT_CGM)
INCLUDES = -I$(top_srcdir)/geom \
           -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom/ACIS \
           -I$(top_srcdir)/geom/virtual \
           -I$(top_srcdir)/geom/facet \
           -I$(top_srcdir)/geom/cholla \
           -I$(srcdir)

noinst_LTLIBRARIES = libcgma_init.la
libcgma_init_la_includedir = $(includedir)
libcgma_init_la_SOURCES = InitCGMA.cpp
libcgma_init_la_include_HEADERS = InitCGMA.hpp


