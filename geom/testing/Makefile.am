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


SUFFIXES += .acis .occ
check_PROGRAMS = $(TESTS) cgm_test
cgm_test_SOURCES = cgm_test.cpp \
                   cgm_test.hpp \
                   TestUtilities.hpp \
                   TestUtilities.cpp \
                   TestConfig.h

%.acis %.occ: cgm_test_script
	cp -f cgm_test_script $@
