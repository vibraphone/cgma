
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
	   $(OCC_INC_FLAG) \
	   -I..\$(OCC_INC_FLAG)

TESTS = makept
if build_ACIS
  TESTS += webcut
endif

check_PROGRAMS = $(TESTS)

webcut_SOURCES = webcut.cpp
makept_SOURCES = makept.cpp
opt_libs=
LDADD = ../geom/virtual/libcubit_virtual.la \
        ../geom/facetbool/libcubit_facetbool.la \
	$(opt_libs) \
        ../geom/Cholla/libcholla.la \
        ../geom/facet/libcubit_facet.la \
        ../geom/libcubit_geom.la \
        ../util/libcubit_util.la 

if build_ACIS
  opt_libs += ../geom/ACIS/libcubit_ACIS.la 
endif
if WITH_OCC
  opt_libs += ../geom/OCC/libcubit_OCC.la $(OCC_LIB_FLAG) $(OCC_LIBS)
endif



