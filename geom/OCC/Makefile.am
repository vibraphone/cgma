# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFAULT_INCLUDES=-I.
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN) $(HAVE_OCC_DEF)
INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           $(OCC_INC_FLAG)

# The name of the library to build
noinst_LTLIBRARIES = libcubit_OCC.la

libcubit_OCC_la_LIBADD = -L$(OCC_LIB_DIR) $(OCC_LIBS)

# The directory where headers will be installed
libcubit_OCC_la_includedir = $(includedir)

# The non-template sources
libcubit_OCC_la_SOURCES = \
    OCCShapeAttributeSet.cpp \
    OCCAttribSet.cpp \
    OCCBody.cpp \
    OCCCoEdge.cpp \
    OCCCoFace.cpp \
    OCCCurve.cpp \
    OCCLoop.cpp \
    OCCLump.cpp \
    OCCModifyEngine.cpp \
    OCCPoint.cpp \
    OCCQueryEngine.cpp \
    OCCShell.cpp \
    OCCSurface.cpp \
    OCCDrawTool.cpp

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_OCC_la_include_HEADERS = \
    OCCShapeAttributeSet.hpp \
    OCCAttribSet.hpp \
    OCCBody.hpp \
    OCCCoEdge.hpp \
    OCCCoFace.hpp \
    OCCCurve.hpp \
    OCCGeometryCreator.hpp \
    OCCLoop.hpp \
    OCCLump.hpp \
    OCCModifyEngine.hpp \
    OCCPoint.hpp \
    OCCQueryEngine.hpp \
    OCCShell.hpp \
    OCCSurface.hpp \
    OCCDrawTool.hpp

