AM_CPPFLAGS += -I$(top_srcdir)/util \
               -I$(top_srcdir)/geom \
               -I$(top_srcdir)/geom/ACIS \
	       -I$(ACIS_DIR)/include \
               -I$(ACIS_DIR)/amendment \
               $(MPI_INCLUDE) \
	       $(HAVE_OCC_DEF)

if USE_MPI
  AM_CPPFLAGS += -DUSE_MPI
endif

# The name of the library to build
if build_parallel
  noinst_LTLIBRARIES = libcubit_parallel.la
else
  noinst_LTLIBRARIES =
endif

# The directory where headers will be installed
libcubit_parallel_la_includedir = $(includedir)

# The non-template sources
libcubit_parallel_la_SOURCES = \
    CGMProcConfig.cpp \
    CGMParallelComm.cpp \
    CGMFileOptions.cpp \
    CABodies.cpp \
    TDParallel.cpp \
    CGMReadParallel.cpp

# The non-template headers
# If any of these do not need to be installed, move them
# to the _SOURCES list above.
libcubit_parallel_la_include_HEADERS = \
    CGMmpi.h \
    CGMmpi_config.h \
    CGMProcConfig.hpp \
    CGMParallelComm.hpp \
    CGMFileOptions.hpp \
    CABodies.hpp \
    TDParallel.hpp \
    CGMReadParallel.hpp

