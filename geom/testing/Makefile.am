BASE_TESTS = AngleCalc.engine \
             GreateGeometry.engine \
             ReadIgesFile.engine \
             GraphicsData.engine

TESTS=
if build_ACIS
  TESTS += $(BASE_TESTS:.engine=.acis)
endif
if build_OCC
  TESTS += $(BASE_TESTS:.engine=.occ)
endif

INCLUDES = -I$(top_srcdir)/util -I$(top_srcdir)/geom -I.
LDADD = ../../util/libcubit_util.la \
        ../libcubit_geom.la

SUFFIXES = .acis .occ
check_PROGRAMS = $(TESTS) cgm_test
cgm_test_SOURCES = cgm_test.cpp \
                   cgm_test.hpp \
                   TestUtilities.hpp \
                   TestUtilities.cpp \
                   TestConfig.h

%.acis %.occ: cgm_test_script
	cp -f cgm_test_script $@
