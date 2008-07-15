AUTOMAKE_OPTIONS = foreign
SUBDIRS = util geom init itaps . test cgm_apps

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = cgm.make

lib_LTLIBRARIES = libcgm.la

EXT_LINK_FLAGS =
EXT_LIBTOOL_FLAGS = 
EXT_LIBS = 
EXT_INCLUDE = 
if build_parallel
  CGM_ll_LIB = geom/parallel/libcubit_parallel.la
else
  CGM_ll_LIB =
endif
if build_ACIS
  CGM_ACIS_LIB = $(CGM_ll_LIB) geom/ACIS/libcubit_ACIS.la 
  EXT_LINK_FLAGS += -L$(ACIS_LIB_DIR)
  EXT_LIBTOOL_FLAGS += -R$(ACIS_LIB_DIR)
  EXT_LIBS += $(ACIS_LIBS)
  EXT_INCLUDE += $(ACIS_INCLUDES) $(ACIS_DEFINES)
else
  CGM_ACIS_LIB = 
endif
if WITH_OCC
  CGM_OCC_LIB = geom/OCC/libcubit_OCC.la
  EXT_LINK_FLAGS += -L$(OCC_LIB_DIR)
  EXT_LIBTOOL_FLAGS += -R$(OCC_LIB_DIR)
  EXT_LIBS += $(OCC_LIBS)
  EXT_INCLUDE += $(OCC_INC_FLAG)
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
  CORE_CGM_LIBS = $(CUBIT_LTFLAGS) $(CUBIT_LDFLAGS) $(CUBIT_LIBS)
  EXT_LINK_FLAGS += $(CUBIT_LDFLAGS)
  EXT_LIBTOOL_FLAGS += $(CUBIT_LTFLAGS)
endif



libcgm_la_LIBADD = init/libcgma_init.la \
                   itaps/libiGeom.la \
                   $(CORE_CGM_LIBS) \
                   $(EXT_LINK_FLAGS) \
                   $(EXT_LIBS)

libcgm_la_SOURCES = 


# By default, cmg.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
CFG_FILE = $(DESTDIR)$(cfgdir)/cgm.make
install-data-hook:
	echo "CGM_EXT_INCLUDE = $(EXT_INCLUDE)"       >> $(CFG_FILE)
	echo "CGM_EXT_LIBS    = $(EXT_LIBS)"          >> $(CFG_FILE)
	echo "CGM_EXT_LDFLAGS = $(EXT_LINK_FLAGS)"    >> $(CFG_FILE)
	echo "CGM_EXT_LTFLAGS = $(EXT_LIBTOOL_FLAGS)" >> $(CFG_FILE)
	echo "CGM_LIBDIR      = $(libdir)"            >> $(CFG_FILE)
	echo "CGM_INCLUDEDIR  = $(includedir)"        >> $(CFG_FILE)
