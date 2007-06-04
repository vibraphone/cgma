# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) \
       $(LITTLE_ENDIAN) $(ACIS_STEP_TRANSLATOR) $(ACIS_IGES_TRANSLATOR) \
       -DCUBIT_ACIS_VERSION=$(ACIS_VERSION) -DACIS_VERSION=$(ACIS_VERSION) \
       -D$(ACIS_PLATFORM)

LIBS := -L$(ACIS_LIB_DIR) $(ACIS_LIBS) $(LIBS)

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           -I$(top_srcdir)/geom/ACIS \
           -I$(ACIS_DIR)/include \
           $(MPI_INCLUDE)

# The name of the library to build
lib_LTLIBRARIES = libcubit_parallel.la

# The directory where headers will be installed
libcubit_parallel_la_includedir = $(includedir)

# The non-template sources
libcubit_parallel_la_SOURCES = \
    AcisMemFile.cpp \
    CABodies.cpp \
    CGMMemFile.cpp \
    ParallelGeomTool.cpp \
    ProcData.cpp \
    TDParallel.cpp

# The non-template headers
# If any of these do not need to be installed, move them
# to the _SOURCES list above.
libcubit_parallel_la_include_HEADERS = \
    AcisMemFile.hpp \
    CABodies.cpp \
    CGMMemFile.hpp \
    ParallelGeomTool.hpp \
    ProcData.hpp \
    TDParallel.hpp

