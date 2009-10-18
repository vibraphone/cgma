AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES = iGeom_protos.h

DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF)

# because "ISO C++ does not support ‘long long"
CXXFLAGS += -Wno-long-long
#CXXFLAGS = -Wall -pipe -pedantic -g -DDEBUG -Wno-long-long -DUSE_MPI

if USE_BABEL
  SIDL_DIR = SIDL
else
  SIDL_DIR = 
endif 

TESTS = 
if build_ACIS
  TESTS += testgeom
else
if WITH_CUBIT
  TESTS += testgeom
endif
endif
if build_OCC
  TESTS += testgeom_occ
endif
if build_parallel
  TESTS += partest
endif


SUBDIRS = . $(SIDL_DIR)

lib_LTLIBRARIES = libiGeom.la

libiGeom_la_includedir = $(includedir)

libiGeom_la_LIBADD = $(top_builddir)/libcgm.la

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
	   -I$(top_srcdir)/init \
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
	iGeom_FCDefs.h \
	iGeom_protos.h \
	iBase_f.h \
	iGeom_f.h

check_PROGRAMS = chaman $(TESTS)

testgeom_SOURCES = testgeom.cc
testgeom_LDADD = libiGeom.la
testgeom_CPPFLAGS = $(CPPFLAGS) -DSRCDIR=$(srcdir)
testgeom_LDFLAGS = $(CGM_EXT_LTFLAGS) $(CGM_EXT_LDFLAGS)
testgeom_occ_SOURCES = $(testgeom_SOURCES)
testgeom_occ_LDADD = $(testgeom_LDADD)
testgeom_occ_CPPFLAGS = $(testgeom_CPPFLAGS) -DFORCE_OCC
testgeom_occ_LDFLAGS = $(testgeom_LDFLAGS)

chaman_SOURCES = chaman.cc
chaman_LDADD = libiGeom.la
chaman_CPPFLAGS = $(CPPFLAGS) -DSRCDIR=$(srcdir)
chaman_LDFLAGS = $(CGM_EXT_LTFLAGS) $(CGM_EXT_LDFLAGS)

if build_parallel
  partest_SOURCES = partest.cpp
  partest_LDADD = $(testgeom_LDADD)
  partest_CPPFLAGS = $(testgeom_CPPFLAGS) -DFORCE_OCC
  partest_LDFLAGS = $(testgeom_LDFLAGS)
  INCLUDES += -I$(top_srcdir)/geom/parallel
endif


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
CFG_FILE = $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
install-data-hook:
	echo "IGEOM_CXX_LDFLAGS = -L${libdir}" >> $(CFG_FILE)
	echo "IGEOM_CXX_LTFLAGS = -R${libdir}" >> $(CFG_FILE)
	echo "IGEOM_INCLUDES = -I${includedir}" >> $(CFG_FILE)

EXTRA_DIST = testgeom.sat size.sat

CLEANFILES = iGeom-Defs.inc testout.sat mmgr.log
DISTCLEANFILES = iGeom_FCDefs.h

FPFX=iGeom
SEDEXPR = 's/^[[:space:]]*void[[:space:]][[:space:]]*$(FPFX)_\([a-z][_a-zA-Z0-9]*\)[[:space:]]*(.*$$/$(FPFX)_\1/p'
iGeom_protos.h: iGeom.h
	echo '#include "iGeom_FCDefs.h"' > $@
	echo '#ifdef IGEOM_FC_FUNC_' >> $@
	echo >> $@
	for func in `$(SED) -n $(SEDEXPR) $<`; do \
	  lower=`echo $$func | tr '[:upper:]' '[:lower:]'`; \
          upper=`echo $$func | tr '[:lower:]' '[:upper:]'`; \
	  echo "#define $$func IGEOM_FC_FUNC_( $$lower, $$upper )" >> $@; \
        done
	echo >> $@
	echo "#endif" >> $@

iGeom_FCDefs.h: iBase_FCDefs.h
	cd .. && ./config.status itaps/iGeom_FCDefs.h

