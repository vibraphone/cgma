AUTOMAKE_OPTIONS = foreign
SUBDIRS = util geom init itaps test cgm_apps

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = cgm.make

# By default, cmg.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo 'CGM_LIBS_LDFLAGS=$$(CGM_LIBS_LDFLAGS_INSTALLED)' >> $(DESTDIR)$(cfgdir)/cgm.make
	echo 'CGMA_INIT_LDFLAGS=$$(CGMA_INIT_LDFLAGS_INSTALLED)' >> $(DESTDIR)$(cfgdir)/cgm.make
	echo 'CGM_INCLUDES=$$(CGM_INC_CPPFLAGS_INSTALLED)' >> $(DESTDIR)$(cfgdir)/cgm.make
	echo "CGM_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/cgm.make
	echo "CGM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/cgm.make
