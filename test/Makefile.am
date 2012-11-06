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

TESTS = 
if build_ACIS
  TESTS += webcut hollow_acis brick_acis merge_acis AngleCalc_acis CreateGeometry_acis GraphicsData_acis
else
if WITH_CUBIT
  TESTS += webcut merge_acis
else
  TESTS += brick_facet merge_facet init sheet
endif
endif
if build_OCC
  TESTS += init sheet attribute_to_file loft offset_curves point_project imprint_bug modify subtract makept test_occ brick_occ merge_occ r_w operation section AngleCalc_occ  CreateGeometry_occ GraphicsData_occ
endif

check_PROGRAMS = $(TESTS)

if BUILD_CGM
  LDADD = ../libcgm.la ../geom/libcubit_geom.la ../init/libcgma_init.la ../util/libcubit_util.la libcgm_test.la ../geom/virtual/libcubit_virtual.la ../geom/facet/libcubit_facet.la $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
else
  LDADD = -lcubit_geom -lcubit_util -lcubiti19 $(top_builddir)/init/libcgma_init.la $(top_builddir)/geom/libcubit_geom.la $(top_builddir)/util/libcubit_util.la $(CGM_EXT_LIBS) $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)
endif

noinst_LTLIBRARIES = libcgm_test.la
libcgm_test_la_SOURCES = TestUtilities.hpp \
                         TestUtilities.cpp

webcut_SOURCES = webcut.cpp
hollow_acis_SOURCES = hollow_acis.cpp
makept_SOURCES = makept.cpp
offset_curves_SOURCES = offset_curves.cpp
loft_SOURCES = loft.cpp
modify_SOURCES = modify.cpp
imprint_bug_SOURCES = imprint_bug.cpp
subtract_SOURCES = subtract.cpp
section_SOURCES = section.cpp
sheet_SOURCES  = sheet.cpp
r_w_SOURCES    = r_w.cpp
test_occ_SOURCES = test_occ.cpp
point_project_SOURCES = point_project.cpp
operation_SOURCES = operation.cpp
init_SOURCES = init.cpp
attribute_to_file_SOURCES = attribute_to_file.cpp
attribute_to_file_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OCC

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

AngleCalc_acis_SOURCES = AngleCalc.cpp
AngleCalc_acis_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_ACIS
AngleCalc_occ_SOURCES = AngleCalc.cpp
AngleCalc_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OCC

CreateGeometry_acis_SOURCES = CreateGeometry.cpp
CreateGeometry_acis_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_ACIS
CreateGeometry_occ_SOURCES = CreateGeometry.cpp
CreateGeometry_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OCC

GraphicsData_acis_SOURCES = GraphicsData.cpp
GraphicsData_acis_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_ACIS
GraphicsData_occ_SOURCES =  GraphicsData.cpp
GraphicsData_occ_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS) -DTEST_OCC

# Files that are used as input to tests, and therefore
# must be included in tarball for tests to run.
EXTRA_DIST = 62_shaver1.brep \
             66_shaver3.brep \
             72_shaver6.brep \
             LeverArm.brep \
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
             webcut.brep     \
             Line_1.brep    \
             Line_2.brep    \
             stitch.name_occ \
             box-w-hole.brep \
             Solid_6.brep  \
             Solid_7.brep  

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
 

