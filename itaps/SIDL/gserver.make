include babel.make

TARGET_LIB = libiGeomserver.la
LTLIBS = ../../libiGeom.la
C_SRCS = $(IORSRCS)
CXX_SRCS = $(IMPLSRCS) $(SKELSRCS) $(STUBSRCS)
LIBLINK = $(CXXLINK)
bdir = gserver

include ../common.make
