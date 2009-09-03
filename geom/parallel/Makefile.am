# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign


# Override default defines with the ones we want from the configure script
if build_OCC
#  DEFS = $(TEMPLATE_DEFS_INCLUDED) $(HAVE_OCC_DEF)
   DEFS += $(HAVE_OCC_DEF)
else
   DEFS += $(TEMPLATE_DEFS_INCLUDED) \
         $(LITTLE_ENDIAN) $(ACIS_STEP_TRANSLATOR) $(ACIS_IGES_TRANSLATOR) \
         -DCUBIT_ACIS_VERSION=$(ACIS_VERSION) -DACIS_VERSION=$(ACIS_VERSION) \
         -D$(ACIS_PLATFORM)
#   DEFS += $(ACIS_STEP_TRANSLATOR) $(ACIS_IGES_TRANSLATOR) \
         -DCUBIT_ACIS_VERSION=$(ACIS_VERSION) -DACIS_VERSION=$(ACIS_VERSION) \
         -D$(ACIS_PLATFORM)
endif

# because "ISO C++ does not support ‘long long"
CXXFLAGS += -Wno-long-long

INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           -I$(top_srcdir)/geom/ACIS \
	   -I$(ACIS_DIR)/include \
           $(MPI_INCLUDE)

# The name of the library to build
noinst_LTLIBRARIES = 
if BUILD_CGM
  noinst_LTLIBRARIES += libcubit_parallel.la
endif

# The directory where headers will be installed
libcubit_parallel_la_includedir = $(includedir)

# The non-template sources
libcubit_parallel_la_SOURCES = \
    CGMProcConfig.cpp \
    CGMParallelComm.cpp \
    FileOptions.cpp \
    ParallelGeomTool.cpp

# The non-template headers
# If any of these do not need to be installed, move them
# to the _SOURCES list above.
libcubit_parallel_la_include_HEADERS = \
    CGMProcConfig.hpp \
    CGMParallelComm.hpp \
    FileOptions.hpp \
    ParallelGeomTool.hpp

