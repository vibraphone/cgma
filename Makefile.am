AUTOMAKE_OPTIONS = foreign
SUBDIRS = util geom 

if build_ACIS
  ACIS_EXAMPLE_TEST=cd cgm_apps/examples/driverc++; make check-TESTS
endif

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = cgm.make

EXTRA_DIST = \
	cgm_apps/examples/README \
	cgm_apps/examples/facetdriver/Makefile \
	cgm_apps/examples/facetdriver/facets.cpp \
	cgm_apps/examples/facetdriver/README \
	cgm_apps/examples/driverc++/Makefile \
	cgm_apps/examples/driverc++/mergechk.cpp \
	cgm_apps/examples/driverc++/merge.sat

example_test:
	cd cgm_apps/examples/facetdriver; make check-TESTS
	$(ACIS_EXAMPLE_TEST)

check: example_test

# By default, cmg.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "CGM_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/cgm.make
	echo "CGM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/cgm.make
