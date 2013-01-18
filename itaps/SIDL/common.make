

LIBS = $(LTLIBS) -L$(BABEL_DIR)/lib -lsidl

INCLUDES = -I$(BABEL_DIR)/include \
	   -I$(builddir) -I$(srcdir) \
           -I$(igeom_srcdir) \
           -I$(igeom_builddir) \
           -I$(sidl_srcdir) \
           -I$(top_srcdir) \
           -I$(top_builddir)

igeom_srcdir = $(sidl_srcdir)/..
igeom_builddir = ../..
sidl_srcdir = ../.
sidl_builddir = ..
srcdir = .././$(bdir)
top_srcdir = ../../..

pkgdatadir = $(datadir)/cgma
pkglibdir = $(libdir)/cgma
pkgincludedir = $(includedir)/cgma
top_builddir = ../../..
builddir = .
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
INSTALL = /usr/bin/install -c
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
target_triplet = x86_64-unknown-linux-gnu
subdir = itaps/SIDL/$(bdir)
libLTLIBRARIES_INSTALL = $(INSTALL)
LTLIBRARIES = $(lib_LTLIBRARIES)

CXXCOMPILE = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
LTCXXCOMPILE = $(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CXXFLAGS) $(CXXFLAGS)
CXXLD = $(CXX)
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --tag=CC --mode=compile $(CC) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(LIBTOOL) --tag=CC --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@

ACLOCAL = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/missing --run aclocal-1.11
AMDEP_FALSE = #
AMDEP_TRUE = 
AMTAR = $${TAR-tar}
AR = ar
AUTOCONF = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/missing --run autoconf
AUTOHEADER = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/missing --run autoheader
AUTOMAKE = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/missing --run automake-1.11
AWK = gawk
BABEL_DIR = 
BABEL = $(BABEL_DIR)/bin/babel
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS =  -Wall -pipe -pedantic -g
CPP = gcc -E
CPPFLAGS = 
CXX = g++
CXXCPP = g++ -E
CXXDEPMODE = depmode=gcc3
CXXFLAGS =  -Wall -pipe -pedantic -g
CYGPATH_W = echo

# Some variables
DEFS = 
DEPDIR = .deps
DISTCHECK_CONFIGURE_FLAGS =  --with-occ="/mnt/disk2b/jhu/build653"
ECHO = @ECHO@
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
GREP = /bin/grep
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LDFLAGS = 
LIBOBJS = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LN_S = ln -s
MAKEINFO = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/missing --run makeinfo
OBJEXT = o
PACKAGE = cgma
PACKAGE_BUGREPORT = 
PACKAGE_NAME = CGMA
PACKAGE_STRING = CGMA 13.1
PACKAGE_TARNAME = cgma
PACKAGE_VERSION = 13.1
PARALLEL_FALSE = @PARALLEL_FALSE@
PARALLEL_HDF5_FALSE = @PARALLEL_HDF5_FALSE@
PARALLEL_HDF5_TRUE = @PARALLEL_HDF5_TRUE@
PARALLEL_TRUE = @PARALLEL_TRUE@
PATH_SEPARATOR = :
RANLIB = ranlib
SET_MAKE = 
SHELL = /bin/bash
STRIP = strip
VERSION = 13.1
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /mnt/disk2b/jhu/merge-cubit14.0/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = /bin/mkdir -p
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /usr/local
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
sysconfdir = ${prefix}/etc
target = x86_64-unknown-linux-gnu
target_alias = 
target_cpu = x86_64
target_os = linux-gnu
target_vendor = unknown

# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# The list of source files, and any header files that do not need to be installed
SOURCES = $(C_SRCS) $(CXX_SRCS)
OBJECTS = $(C_SRCS:.c=.lo) $(CXX_SRCS:.cc=.lo)
HEADERS = $(IMPLHDRS) $(IORHDRS) $(STUBHDRS)

all: $(TARGET_LIB)

install: install-lib install-headers

uninstall: uninstall-lib uninstall-headers

$(TARGET_LIB): $(OBJECTS) $(LTLIBS)
	$(LIBLINK) -rpath $(libdir) $(LDFLAGS) $(OBJECTS) $(LIBADD) $(LIBS)


install-lib: $(TARGET_LIB)
	test -d "$(libdir)" || mkdir "$(libdir)"
	$(LIBTOOL) --mode=install $(libLTLIBRARIES_INSTALL) $(INSTALL_STRIP_FLAG) $< '$(DESTDIR)$(libdir)/$<'

uninstall-lib:
	$(LIBTOOL) --mode=uninstall rm -f "$(DESTDIR)$(libdir)/$(TARGET_LIB)"

install-headers: $(HEADERS)
	test -d "$(includedir)" || mkdir "$(includedir)"
	@list='$(HEADERS)'; for p in $$list; do \
	  $(ECHO) $(INSTALL_HEADER) $$p $(DESTDIR)$(includedir)/$$p ; \
	  $(INSTALL_HEADER) "$$p" "$(DESTDIR)$(includedir)/$$p" ; \
	done

uninstall-headers:
	@list='$(HEADERS)'; for p in $$list; do \
	  $(ECHO) "rm -f $(DESTDIR)$(includedir)/$$p" ; \
	  rm -f "$(DESTDIR)$(includedir)/$$p" ; \
	done

.cc.lo:
	$(LTCXXCOMPILE) -c -o $@ $<

.c.lo:
	$(LTCOMPILE) -c -o $@ $<

.SUFFIXES: .lo .cc .c

.PHONEY: all install uninstall install-lib uninstall-lib install-headers uninstall-headers

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:


# dependencies
#
# babel generates a babel.make.dependencies, but at least for 
# 0.10.10, everything it lists depends either directly or 
# indirectly on the SIDL file.  As we're using the 'babel.make'
# as a timestamp, just depend on that to achieve nearly the 
# same effect.
$(OBJECTS) : babel.make

# server implementation files also depend on the C interface
IGEOM_C_INC = $(igeom_srcdir)/iBase.h \
              $(igeom_srcdir)/iGeom.h \
	      $(igeom_srcdir)/iGeom_protos.h

IMPLOBJS = $(IMPLSRCS:.cc=.lo)
$(IMPLOBJS) : $(IGEOM_C_INC)
