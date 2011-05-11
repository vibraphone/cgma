BASE_TESTS = AngleCalc.engine \
             CreateGeometry.engine \
             GraphicsData.engine
#             ReadIgesFile.engine 

TESTS=
if build_ACIS
  TESTS += $(BASE_TESTS:.engine=.acis)
  AM_CPPFLAGS += -DCGM_ACIS -I$(srcdir)/../ACIS
endif
if build_OCC
  TESTS += $(BASE_TESTS:.engine=.occ)
  AM_CPPFLAGS += -DCGM_OCC -I$(srcdir)/../OCC $(OCC_INC_FLAG)
endif

AM_CPPFLAGS += -I$(top_srcdir)/util -I$(top_srcdir)/geom -I.
LDADD = libcgm_test.la ../../libcgm.la  $(CGM_EXT_LDFLAGS) $(CGM_EXT_LTFLAGS)

check_PROGRAMS = $(BASE_TESTS:.engine=)
check_SCRIPT = $(TESTS)

noinst_LTLIBRARIES = libcgm_test.la
libcgm_test_la_SOURCES = TestConfig.h \
                      cgm_test.cpp \
                      cgm_test.hpp \
                      TestUtilities.hpp \
                      TestUtilities.cpp

AngleCalc_SOURCES = main.cpp AngleCalc.cpp
AngleCalc_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) -DTEST=AngleCalc
CreateGeometry_SOURCES = main.cpp CreateGeometry.cpp
CreateGeometry_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) -DTEST=CreateGeometry
#ReadIgesFile_SOURCES = main.cpp ReadIgesFile.cpp
#ReadIgesFile_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) -DTEST=ReadIgesFile
GraphicsData_SOURCES = main.cpp GraphicsData.cpp
GraphicsData_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) -DTEST=GraphicsData

$(TESTS): Makefile
	$(AM_V_GEN) test=`expr $@ : '\(.*\)\.'` ;\
	engine=`expr $@ : '.*\.\(.*\)'` ;\
	echo "#!/bin/sh" >$@ ;\
	echo "./$$test -E $$engine" >>$@ ;\
	chmod a+x $@

MOSTLYCLEANFILES = mmgr.log
CLEANFILES = $(TESTS)
