
# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF) $(CUBIT_CGM)
INCLUDES = -I$(top_srcdir)/geom \
           -I$(top_srcdir)/util \
           -I$(top_srcdir)/init \
           -I$(top_srcdir)/geom/ACIS \
           -I$(top_srcdir)/geom/virtual \
           -I$(top_srcdir)/geom/facet \
           -I$(top_srcdir)/geom/Cholla \
	   -I$(top_srcdir)/geom/OCC \
	   -I$(srcdir) \
	   $(OCC_INC_FLAG)

TESTS = sheet
if build_ACIS
  TESTS += webcut hollow_acis
else
if WITH_CUBIT
  TESTS += webcut
endif
endif
if build_OCC
  TESTS += modify makept r_w test_occ
endif

AM_CPPFLAGS =
AM_CPPFLAGS += -DSRCDIR=@srcdir@ $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF)
check_PROGRAMS = $(TESTS)

LINK_FLAGS = $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)

webcut_SOURCES = webcut.cpp
hollow_acis_SOURCES = hollow_acis.cpp
makept_SOURCES = makept.cpp
modify_SOURCES = modify.cpp
sheet_SOURCES  = sheet.cpp
r_w_SOURCES    = r_w.cpp
test_occ_SOURCES = test_occ.cpp
LDADD = ../libcgm.la
webcut_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
hollow_acis_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
makept_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
modify_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
sheet_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
r_w_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)
test_occ_LDFLAGS = $(LDFLAGS) $(LINK_FLAGS)

EXTRA_DIST = 62_shaver1.brep \
             66_shaver3.brep \
             72_shaver6.brep \
             LeverArm.brep \
             aftersub.occ_ref \
             beforesub.occ_ref \
             beforesub2.occ_ref \
             unite1.occ\
             unite2.occ\
             unite3.occ\
             unite4.occ\
             unite5.occ\
             unite6.occ\
             huge.sat \
             model1.sat \
             model2.sat \
             model3.sat \
             model4.sat \
             model5.sat \
             model6.sat \
             model7.sat \
             hollow.sat \
             stitch.occ_name \
             webcut_with_brick.sat \
             webcut_with_curve_loop.sat \
             webcut_with_extended_surf.sat \
             webcut_with_sheet.sat \
             webcut_with_sweep_along_curve.sat \
             webcut_with_sweep_surfaces.sat \
             webcut_with_sweep_surfaces_along_curve.sat \
             webcut_with_sweep_surfaces_rotated.sat

