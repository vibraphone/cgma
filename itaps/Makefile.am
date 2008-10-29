AUTOMAKE_OPTIONS = foreign

DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_ACIS_DEF)

if USE_BABEL
  SIDL_DIR = SIDL
else
  SIDL_DIR = 
endif 

SUBDIRS = . $(SIDL_DIR)

lib_LTLIBRARIES = libiGeom.la

libiGeom_la_includedir = $(includedir)

libiGeom_la_LIBADD = ../libcgm.la

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
	   -I$(top_srcdir)/init \
           -I$(top_srcdir)/geom/ACIS \
           -I.

libiGeom_la_SOURCES = \
        CGMAIterator.hpp \
        CATag.hpp \
	CATag.cpp \
	iGeom_CGMA.cc \
	iGeomError.cc \
	iGeomError.h 

libiGeom_la_include_HEADERS = \
	iBase.h \
	iGeom.h \
	iBase_FCDefs.h \
	iGeom_protos.h \
	iBase_f.h \
	iGeom_f.h

TESTS = testgeom

check_PROGRAMS = $(TESTS)

testgeom_SOURCES = testgeom.cc
testgeom_LDADD = libiGeom.la
testgeom_CPPFLAGS = $(CPPFLAGS) -DSRCDIR=$(srcdir)
testgeom_LDFLAGS = $(CGM_LIBS_LTFLAGS)

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
	echo "include $(libdir)/cgm.make"        >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_LDFLAGS  = $$(CGM_LDFLAGS) $$(IGEOM_CXX_LDFLAGS)'  >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_LTFLAGS  = $$(CGM_LTFLAGS)'  >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_LIBS     = $$(CGM_LIBS) $$(IGEOM_CXX_LIBS)'     >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_DEFINES  = $$(CGM_DEFINES)'  >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_INCLUDES = $$(CGM_INCLUDES)' >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo 'IGEOM_CPPFLAGS = $$(CGM_CPPFLAGS)' >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc

EXTRA_DIST = testgeom.sat

CLEANFILES = iGeom-Defs.inc
