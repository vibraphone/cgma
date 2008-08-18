AUTOMAKE_OPTIONS = foreign

DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_ACIS_DEF)

if USE_BABEL
  SIDL_DIR = SIDL
else
  SIDL_DIR = 
endif 

SUBDIRS = . $(SIDL_DIR)

noinst_LTLIBRARIES = libiGeom.la

if build_ACIS
  TEST_ACISPATH = -R$(ACIS_LIB_DIR)
else
  TEST_ACISPATH =
endif

libiGeom_la_includedir = $(includedir)

libiGeom_la_libadd = $(top_builddir)/libcgm.la

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
	iGeom_protos.h 

TESTS = testgeom

check_PROGRAMS = $(TESTS)

testgeom_SOURCES = testgeom.cc
testgeom_LDADD = ../libcgm.la
testgeom_CPPFLAGS = $(CPPFLAGS) -DSRCDIR=$(srcdir)
testgeom_LDFLAGS = $(TEST_ACISPATH)

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
iGeom-Defs.inc: Makefile.am
	@echo "include $(cfgdir)/cgm.make"         > $@
	@echo 'iGeom_LDFLAGS  = $$(CGM_LDFLAGS)'  >> $@
	@echo 'iGeom_LTFLAGS  = $$(CGM_LTFLAGS)'  >> $@
	@echo 'iGeom_LIBS     = $$(CGM_LIBS)'     >> $@
	@echo 'iGeom_DEFINES  = $$(CGM_DEFINES)'  >> $@
	@echo 'iGeom_INCLUDES = $$(CGM_INCLUDES)' >> $@
	@echo 'iGeom_CPPFLAGS = $$(CGM_CPPFLAGS)' >> $@

CLEANFILES = iGeom-Defs.inc
