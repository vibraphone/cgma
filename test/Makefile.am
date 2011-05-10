AM_CPPFLAGS += $(HAVE_ACIS_DEF) $(HAVE_OCC_DEF) $(CUBIT_CGM) \
           -DSRCDIR=@srcdir@ \
           -I$(top_srcdir)/geom \
           -I$(top_srcdir)/util \
           -I$(top_srcdir)/init \
           -I$(top_srcdir)/geom/ACIS \
           -I$(top_srcdir)/geom/virtual \
           -I$(top_srcdir)/geom/facet \
           -I$(top_srcdir)/geom/Cholla \
	   -I$(top_srcdir)/geom/OCC \
	   -I$(srcdir) \
	   $(OCC_INC_FLAG)

TESTS = init sheet
if build_ACIS
  TESTS += webcut hollow_acis brick_acis merge_acis
else
if WITH_CUBIT
  TESTS += webcut brick_acis merge_acis
else
  TESTS += brick_facet merge_facet
endif
endif
if build_OCC
  TESTS += modify subtract makept r_w test_occ operation brick_occ merge_occ
endif

check_PROGRAMS = $(TESTS)

LDADD = ../libcgm.la
AM_LDFLAGS += $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)

webcut_SOURCES = webcut.cpp
hollow_acis_SOURCES = hollow_acis.cpp
makept_SOURCES = makept.cpp
modify_SOURCES = modify.cpp
subtract_SOURCES = subtract.cpp
sheet_SOURCES  = sheet.cpp
r_w_SOURCES    = r_w.cpp
test_occ_SOURCES = test_occ.cpp
operation_SOURCES = operation.cpp
init_SOURCES = init.cpp

brick_sources = brick.cpp
brick_acis_SOURCES = $(brick_sources)
brick_acis_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) '-DTEST_ENGINE="ACIS"'
brick_occ_SOURCES = $(brick_sources)
brick_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) '-DTEST_ENGINE="OCC"'
brick_facet_SOURCES = $(brick_sources)
brick_facet_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) '-DTEST_ENGINE="FACET"'

merge_sources = merge_test.cpp
merge_acis_SOURCES = $(merge_sources)
merge_acis_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_ACIS
merge_occ_SOURCES = $(merge_sources)
merge_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OCC
merge_facet_SOURCES = $(merge_sources)
merge_facet_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_FACET

# Files that are used as input to tests, and therefore
# must be included in tarball for tests to run.
EXTRA_DIST = 62_shaver1.brep \
             66_shaver3.brep \
             72_shaver6.brep \
             LeverArm.brep \
             aftersub.occ_ref \
             beforesub.occ_ref \
             beforesub2.occ_ref \
             unite1.occ\
             huge.sat \
             model1.sat \
             model2.sat \
             model3.sat \
             model4.sat \
             model5.sat \
             model6.sat \
             model7.sat \
             hollow.sat \
	     unite4.occ \
             proe.stp   \
             diffuser.iges   \
             ex3.iges   \
             Cylinder_1.brep   \
             Cylinder_2.brep   \
             stitch.name_occ 

# Files that are created by tests, and that therefore 
#  - should be removed during 'make clean'
#  - should not be included in tarball
MOSTLYCLEANFILES = \
	webcut_with_brick.sat \
	webcut_with_cylinder.sat \
	webcut_with_sheet.sat \
	webcut_with_sweep_rotational.sat \
	webcut_with_sweep_curves.sat \
	webcut_with_sweep_along_curve.sat \
	webcut_with_curve_loop.sat \
	webcut_with_extended_surf.sat \
	webcut_with_sweep_surfaces_rotated.sat \
	webcut_with_sweep_surfaces.sat \
	webcut_with_sweep_surfaces_along_curve.sat \
	webcut_with_sweep_surfaces_perp.sat \
	point.occ \
	stitch.occ \
	beforesub.occ \
	beforesub2.occ \
	aftersub.occ \
	unite2.occ \
	unite3.occ \
	unite5.occ \
	unite6.occ \
	proeexport.iges \
	proe.occ \
	diffuser.occ \
	proeexport.step \
	ex3export.step \
	ex3.occ \
	ex3export.iges \
        prism.occ \
        pyramid.occ \
        fillet.occ \
	mmgr.log
 

