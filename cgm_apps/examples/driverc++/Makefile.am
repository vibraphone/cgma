TESTS = mergechk_merge mergechk_huge
EXTRA_DIST = README merge.sat huge.sat $(TESTS)

check_PROGRAMS = mergechk
CLEANFILES = mmgr.log

AM_CPPFLAGS += -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/ACIS

AM_LDFLAGS = -L$(ACIS_LIB_DIR) -R$(ACIS_LIB_DIR)

mergechk_SOURCES = mergechk.cpp

LDADD = ../../../libcgm.la

