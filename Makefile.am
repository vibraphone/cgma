AUTOMAKE_OPTIONS = foreign
SUBDIRS = compat util geom init . itaps test cgm_apps
# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = cgm.make

lib_LTLIBRARIES = libcgm.la

if build_parallel
  CGM_ll_LIB = geom/parallel/libcubit_parallel.la
else
  CGM_ll_LIB =
endif
if build_ACIS
  CGM_ACIS_LIB = $(CGM_ll_LIB) geom/ACIS/libcubit_ACIS.la 
  
else
  CGM_ACIS_LIB = 
endif
if build_OCC
  CGM_OCC_LIB = geom/OCC/libcubit_OCC.la
else
  CGM_OCC_LIB =
endif

if BUILD_CGM
  CORE_CGM_LIBS = $(CGM_ACIS_LIB) \
                  $(CGM_OCC_LIB) \
                  $(CGM_ll_LIB) \
                  geom/virtual/libcubit_virtual.la \
                  geom/facet/libcubit_facet.la \
		  geom/Cholla/libcholla.la \
                  geom/facetbool/libcubit_facetbool.la \
                  geom/libcubit_geom.la \
                  util/libcubit_util.la
else
  CORE_CGM_LIBS = 
endif



libcgm_la_LIBADD = init/libcgma_init.la $(CORE_CGM_LIBS) $(CGM_EXT_LIBS) $(CUBIT_FILE)

libcgm_la_SOURCES = 

# By default, cmg.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
CFG_FILE = $(DESTDIR)$(cfgdir)/cgm.make
install-data-hook:
	$(AM_V_at)echo "CGM_INT_INCLUDE = -I$(includedir)"  >> $(CFG_FILE)
	$(AM_V_at)echo "CGM_INT_LDFLAGS = -L$(libdir)"      >> $(CFG_FILE)
	$(AM_V_at)echo "CGM_INT_LTFLAGS = -R$(libdir)"      >> $(CFG_FILE)
