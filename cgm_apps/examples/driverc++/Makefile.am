TESTS = mergechk_merge mergechk_huge
EXTRA_DIST = README merge.sat huge.sat $(TESTS)

check_PROGRAMS = mergechk

CPPFLAGS += -I$(top_srcdir)/util \
            -I$(top_srcdir)/geom \
            -I$(top_srcdir)/geom/ACIS

LDFLAGS += -L$(ACIS_LIB_DIR) -R$(ACIS_LIB_DIR)

mergechk_SOURCES = mergechk.cpp

LDADD = ../../../geom/libcubit_geom.la \
        ../../../geom/ACIS/libcubit_ACIS.la

