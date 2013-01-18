include babel.make

TARGET_LIB = libiGeomCclient.la
LTLIBS = ../gserver/libiGeomserver.la
C_SRCS = $(STUBSRCS)
CXX_SRCS =
LIBLINK = $(LINK)
bdir = Cclient

include ../common.make
