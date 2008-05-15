AUTOMAKE_OPTIONS = foreign

DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_ACIS_DEF)

if USE_BABEL
  SIDL_DIR = SIDL
else
  SIDL_DIR = 
endif 

SUBDIRS = . $(SIDL_DIR)

lib_LTLIBRARIES = libiGeom.la

if build_ACIS
  TEST_ACISLIBS = $(top_builddir)/geom/ACIS/libcubit_ACIS.la
else
  TEST_ACISLIBS =
endif

if BUILD_CGM
  TEST_XLIBS = $(TEST_ACISLIBS) \
               $(top_builddir)/geom/facet/libcubit_facet.la \
               $(top_builddir)/geom/facetbool/libcubit_facetbool.la \
               $(top_builddir)/geom/Cholla/libcholla.la \
               $(top_builddir)/geom/virtual/libcubit_virtual.la \
               $(top_builddir)/geom/libcubit_geom.la \
               $(top_builddir)/util/libcubit_util.la
else
  TEST_XLDFLAGS = $(CUBIT_LDFLAGS) $(CUBIT_LTFLAGS)
  TEST_XLIBS = $(CUBIT_LIBS)
endif

if build_ACIS
  ACIS_CGM_INCLUDES = -I$(top_srcdir)/geom/ACIS
else
  ACIS_CGM_INCLUDES =
endif

libiGeom_la_includedir = $(includedir)

libiGeom_la_libadd = $(top_builddir)/init/libcgma_init.la $(TEST_XLIBS)

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
	   -I$(top_srcdir)/init \
           -I.

libiGeom_la_SOURCES = \
        CATag.hpp \
	CATag.cpp \
	iGeom_CGMA.cc \
	iGeomError.cc \
	iGeomError.h 

libiGeom_la_include_HEADERS = \
	iBase.h \
	iGeom.h \
	iBase_FCDefs.h \
	iGeom_protos.h 

TESTS = testgeom

check_PROGRAMS = $(TESTS)

testgeom_SOURCES = testgeom.cc
testgeom_LDADD = $(TEST_XLDFLAGS) libiGeom.la  ../init/libcgma_init.la $(TEST_XLIBS)
testgeom_CPPFLAGS = $(CPPFLAGS) -DSRCDIR=$(srcdir)


# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iGeom-Defs.inc

# By default, iGeom-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "iGeom_SRCDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "iGeom_BLDDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "iGeom_DIR=${cfgdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "iGeom_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "CGM_LDFLAGS=-L${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "CGM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc

