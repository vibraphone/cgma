BUILT_SOURCES = iGeom_protos.h

if USE_BABEL
  SIDL_DIR = SIDL
else
  SIDL_DIR = 
endif 

TESTS = test_edge_orient
if build_ACIS
  TESTS += testgeom 
else
if WITH_CUBIT
  TESTS += testgeom 
endif
endif
if build_OCC
  TESTS += testgeom_occ copy_attrib
endif


SUBDIRS = . $(SIDL_DIR)

lib_LTLIBRARIES = libiGeom.la

libiGeom_la_includedir = $(includedir)

if BUILD_CGM
  libiGeom_la_LIBADD = $(top_builddir)/libcgm.la
else
  libiGeom_la_LIBADD = -lcubiti19 -lcubit_geom -lcubit_util $(top_builddir)/init/libcgma_init.la $(top_builddir)/geom/libcubit_geom.la
endif 

AM_CPPFLAGS += $(HAVE_ACIS_DEF) \
               $(HAVE_OCC_DEF) \
               -DSRCDIR=$(srcdir) \
               -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
	       -I$(top_srcdir)/init

if USE_MPI
  AM_CPPFLAGS += -I$(top_srcdir)/geom/parallel -DUSE_MPI
endif

libiGeom_la_SOURCES = \
	CGMAIterator.hpp \
	CATag.hpp \
	CATag.cpp \
	iGeom_CGMA.cc \
	iGeomError.cc \
	iGeomError.h  

libiGeom_la_include_HEADERS = \
	iBase.h \
	iGeom_FCDefs.h \
	iGeom_protos.h \
	iBase_f.h 

if ITAPS_SHIM
  AM_CPPFLAGS += -DITAPS_SHIM
  libiGeom_la_SOURCES += \
    iGeom_private.h
  libiGeom_la_include_HEADERS += \
    CGM_iGeom.h
else
  libiGeom_la_include_HEADERS += \
	iGeom.h \
	iGeom_f.h
endif

check_PROGRAMS = chaman $(TESTS)

if BUILD_CGM
  LDADD = libiGeom.la $(top_builddir)/geom/libcubit_geom.la $(top_builddir)/util/libcubit_util.la $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
else
  LDADD = libiGeom.la $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
endif

AM_LDFLAGS +=  $(CGM_EXT_LTFLAGS) $(CGM_EXT_LDFLAGS)
testgeom_SOURCES = testgeom.cc
testgeom_occ_SOURCES = $(testgeom_SOURCES)
testgeom_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DFORCE_OCC

test_edge_orient_SOURCES = test_edge_orient.cc

chaman_SOURCES = chaman.cc

copy_attrib_SOURCES = copy_attrib.cpp

if USE_MPI
  check_PROGRAMS += partest
  partest_SOURCES = partest.cpp
  partest_CPPFLAGS = $(testgeom_occ_CPPFLAGS)
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
	$(AM_V_at)echo "IGEOM_CXX_LDFLAGS = -L${libdir}" >> $(CFG_FILE)
	$(AM_V_at)echo "IGEOM_CXX_LTFLAGS = -R${libdir}" >> $(CFG_FILE)
	$(AM_V_at)echo "IGEOM_INCLUDEDIR = ${includedir}" >> $(CFG_FILE)
	$(AM_V_at)echo 'IGEOM_INCLUDES = -I$$(IGEOM_INCLUDEDIR)' >> $(CFG_FILE)

EXTRA_DIST = testgeom.sat size.sat

CLEANFILES = iGeom-Defs.inc testout.sat mmgr.log step_import.log test.occ
DISTCLEANFILES = iGeom_FCDefs.h

FPFX=iGeom
SEDEXPR = 's/^[[:space:]]*void[[:space:]][[:space:]]*$(FPFX)_\([a-z][_a-zA-Z0-9]*\)[[:space:]]*(.*$$/$(FPFX)_\1/p'
iGeom_protos.h: iGeom.h
	$(AM_V_GEN)echo '#include "iGeom_FCDefs.h"' > $@ && \
	echo '#ifdef IGEOM_FC_FUNC_' >> $@ && \
	echo >> $@ && \
	for func in `$(SED) -n $(SEDEXPR) $<`; do \
	  lower=`echo $$func | tr '[:upper:]' '[:lower:]'`; \
          upper=`echo $$func | tr '[:lower:]' '[:upper:]'`; \
	  echo "#define $$func IGEOM_FC_FUNC_( $$lower, $$upper )" >> $@; \
        done && \
	echo >> $@ && \
	echo "#endif" >> $@

iGeom_FCDefs.h: iBase_FCDefs.h
	$(AM_V_GEN)cd .. && ./config.status itaps/iGeom_FCDefs.h

