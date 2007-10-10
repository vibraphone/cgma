# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFAULT_INCLUDES=-I.
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)
LIBS := $(OCC_LIB_FLAG) 
INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           -I$(top_srcdir)/geom/Cholla \
           -I$(top_srcdir)/geom/facetbool\
           $(OCC_INC_FLAG)

# The name of the library to build
OCC_TARGET=
if WITH_OCC
  OCC_TARGET += libcubit_OCC.la
endif

lib_LTLIBRARIES = 
if BUILD_CGM
  lib_LTLIBRARIES += $(OCC_TARGET)
endif

# The directory where headers will be installed
libcubit_OCC_la_includedir = $(includedir)

# The non-template sources
libcubit_OCC_la_SOURCES = \
    OCCAttribSet.cpp \
    OCCBody.cpp \
    OCCCoEdge.cpp \
    OCCCurve.cpp \
    OCCLoop.cpp \
    OCCLump.cpp \
    OCCModifyEngine.cpp \
    OCCParamTool.cpp \
    OCCPoint.cpp \
    OCCQueryEngine.cpp \
    OCCShell.cpp \
    OCCSurface.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_OCC_la_include_HEADERS = \
    OCCAttrib.hpp \
    OCCAttribSet.hpp \
    OCCBody.hpp \
    OCCCoEdge.hpp \
    OCCCurve.hpp \
    OCCGeometryCreator.hpp \
    OCCLoop.hpp \
    OCCLump.hpp \
    OCCModifyEngine.hpp \
    OCCParamTool.hpp \
    OCCPoint.hpp \
    OCCQueryEngine.hpp \
    OCCShell.hpp \
    OCCSurface.hpp 

