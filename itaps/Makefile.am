AUTOMAKE_OPTIONS = foreign

DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_ACIS_DEF)

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
               $(top_builddir)/geom/libcubit_geom.la \
               $(top_builddir)/util/libcubit_util.la
  TEST_XLDFLAGS = -R$(ACIS_LIB_DIR)
else
  TEST_XLDFLAGS = -L$(CUBIT_BIND_DIR) -R$(CUBIT_BIN_DIR) 
  TEST_XLIBS = $(CUBIT_LINK)
endif

libiGeom_la_includedir = $(includedir)

if build_ACIS
  ACIS_SRCS = CGMA_ACIS.cc
  ACIS_CGM_INCLUDES = -I$(top_srcdir)/geom/ACIS
else
  ACIS_SRCS = 
  ACIS_CGM_INCLUDES =
endif

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           $(ACIS_CGM_INCLUDES ) \
           -I.

libiGeom_la_SOURCES = \
        CATag.hpp \
	CATag.cpp \
	iGeom_CGMA.cc \
	iGeomError.cc \
	iGeomError.h \
	CGMA_ACIS.h \
	$(ACIS_SRCS)

libiGeom_la_include_HEADERS = \
	iBase.h \
	iGeom.h 

TESTS = testgeom

check_PROGRAMS = $(TESTS)

testgeom_SOURCES = testgeom.cc
testgeom_LDFLAGS = $(TEST_XLDFLAGS)
testgeom_LDADD = libiGeom.la $(TEST_XLIBS)


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
	echo "IGEOM_SRCDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "IGEOM_BLDDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "IGEOM_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "CGM_LDFLAGS=-L${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "CGM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc

