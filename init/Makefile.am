AM_CPPFLAGS += $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF) $(CUBIT_CGM) \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom/ACIS \
               -I$(top_srcdir)/geom/OCC \
               -I$(top_srcdir)/geom/virtual \
               -I$(top_srcdir)/geom/facet \
               -I$(top_srcdir)/geom/cholla \
               -I$(srcdir) \
               $(OCC_INC_FLAG)

noinst_LTLIBRARIES = libcgma_init.la
libcgma_init_la_includedir = $(includedir)
libcgma_init_la_SOURCES = InitCGMA.cpp
libcgma_init_la_include_HEADERS = InitCGMA.hpp


