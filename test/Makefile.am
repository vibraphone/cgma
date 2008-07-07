
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

LINK_FLAGS = 
TESTS = sheet
if build_ACIS
  TESTS += webcut
  LINK_FLAGS += -R$(ACIS_LIB_DIR)
endif
if WITH_OCC
  TESTS += modify makept
  LINK_FLAGS += -R$(OCC_LIB_DIR)
endif

AM_CPPFLAGS =
AM_CPPFLAGS += -DSRCDIR=@srcdir@ $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF)
check_PROGRAMS = $(TESTS)

webcut_SOURCES = webcut.cpp
makept_SOURCES = makept.cpp
modify_SOURCES = modify.cpp
sheet_SOURCES  = sheet.cpp
LDADD = ../libcgm.la
webcut_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
makept_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
modify_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
sheet_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
