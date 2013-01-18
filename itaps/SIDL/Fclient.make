include babel.make

TARGET_LIB = libiGeomFclient.la
LTLIBS = ../gserver/libiGeomserver.la
C_SRCS = $(STUBSRCS)
CXX_SRCS =
LIBLINK = $(LINK)
bdir = Fclient

include ../common.make
