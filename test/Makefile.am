
# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_ACIS_DEF) $(CUBIT_CGM)
INCLUDES = -I$(top_srcdir)/geom \
           -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom/ACIS \
           -I$(top_srcdir)/geom/virtual \
           -I$(top_srcdir)/geom/facet \
           -I$(top_srcdir)/geom/Cholla \
	   -I$(top_srcdir)/geom/OCC \
	   -I$(srcdir) \
	   $(OCC_INC_FLAG)

TESTS=
if build_ACIS
  TESTS += webcut
endif
if WITH_OCC
  TESTS += modify makept
endif

CPPFLAGS += -DSRCDIR=@srcdir@
check_PROGRAMS = $(TESTS)

webcut_SOURCES = webcut.cpp
makept_SOURCES = makept.cpp
modify_SOURCES = modify.cpp
LDADD = ../geom/virtual/libcubit_virtual.la \
        ../geom/facetbool/libcubit_facetbool.la \
        ../geom/Cholla/libcholla.la \
        ../geom/facet/libcubit_facet.la \
        ../geom/libcubit_geom.la \
        ../util/libcubit_util.la 
webcut_LDFLAGS = $(LDFLAGS) -R$(ACIS_LIB_DIR)
webcut_LDADD = ../geom/ACIS/libcubit_ACIS.la  $(LDADD)
makept_LDFLAGS = $(LDFLAGS) -R$(OCC_LIB_DIR)
makept_LDADD = ../geom/OCC/libcubit_OCC.la $(LDADD)
modify_LDFLAGS = $(LDFLAGS) -R$(OCC_LIB_DIR)
modify_LDADD = ../geom/OCC/libcubit_OCC.la $(LDADD)
