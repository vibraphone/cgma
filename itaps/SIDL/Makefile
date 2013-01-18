

MAKE_INC = iGeom-SIDL-Defs.inc
SIDL_FILES = iBase.sidl iGeom.sidl iGeom_SIDL.sidl
SIDL_FILES_PATH = $(srcdir)/iBase.sidl $(srcdir)/iGeom.sidl $(srcdir)/iGeom_SIDL.sidl
# SERVER_FILES = gserver/iBase_Error_Impl.cc \
#                gserver/iBase_Error_Impl.hh \
# 	       gserver/iGeom_SIDL_GeomSidl_Impl.cc \
# 	       gserver/iGeom_SIDL_GeomSidl_Impl.hh \
# 	       gserver/iGeom_Factory_Impl.cc
# SERVER_FILES_PATH = $(srcdir)/gserver/iBase_Error_Impl.cc \
#                     $(srcdir)/gserver/iBase_Error_Impl.hh \
# 	            $(srcdir)/gserver/iGeom_SIDL_GeomSidl_Impl.cc \
# 	            $(srcdir)/gserver/iGeom_SIDL_GeomSidl_Impl.hh 
SERVER_FILES = gserver/iGeom_SIDL_GeomSidl_Impl.cc \
	       gserver/iGeom_SIDL_GeomSidl_Impl.hh \
	       gserver/iGeom_Factory_Impl.cc
SERVER_FILES_PATH = $(srcdir)/gserver/iGeom_SIDL_GeomSidl_Impl.cc \
	            $(srcdir)/gserver/iGeom_SIDL_GeomSidl_Impl.hh 
SOURCE_FILES = $(SIDL_FILES) $(SERVER_FILES) 
DIST_FILES = $(SOURCE_FILES) testgeom.cpp \
             gserver.make Cclient.make Fclient.make
INCLUDES =  -I$(BABEL_DIR)/include -Igserver

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

INSTALL_DATA = ${INSTALL} -m 644
INSTALL = /usr/bin/install -c
BABEL_DIR = 
USE_BABEL = no
BABEL = $(BABEL_DIR)/bin/babel
BABELFLAGS = "--suppress-timestamp"
ECHO = @ECHO@
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
LIBTOOL = $(SHELL) $(top_builddir)/libtool

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
srcdir = .
subdir = itaps/SIDL
sysconfdir = ${prefix}/etc
target = x86_64-unknown-linux-gnu
target_alias = 
target_cpu = x86_64
target_os = linux-gnu
target_vendor = unknown
top_builddir = ../..
top_srcdir = ../..



all: $(MAKE_INC) recursive 
	if test "x$(USE_BABEL)" = "xyes"; then \
	  (cd gserver && $(MAKE) -f ../$(srcdir)/gserver.make $@) && \
	  (cd Cclient && $(MAKE) -f ../$(srcdir)/Cclient.make $@) && \
	  (cd Fclient && $(MAKE) -f ../$(srcdir)/Fclient.make $@) ; \
	fi

clean: 
	rm -rf Cclient Fclient
	rm -rf testgeom testgeom.lo testgeom.o .deps .libs
	if test "x$(srcdir)" != "x."; then \
	  rm -rf gserver; \
	else \
	  for file in gserver/*; do \
	    found=0; \
	    for sfile in $(SERVER_FILES); do \
	      if test "x$$file" = "x$$sfile"; then \
	        found=1; \
              fi ; \
            done ; \
            if test "$$found" -eq 0; then rm -f $$file; fi ; \
	  done ; \
	fi

distclean: clean
	rm -f Makefile iGeom-SIDL-Defs.inc

mostlyclean: 
	for dir in Cclient Fclient gclient; do \
	  rm -f $$dir/*.o $$dir/*.o $$dir/*.la $$dir/*.a $$dir/*.so; \
	  rm -rf $$dir/.libs; \
	done

.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

common.make: $(srcdir)/common.make.in $(top_builddir)/config.status
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

$(MAKE_INC): $(srcdir)/$(MAKE_INC).in $(top_builddir)/config.status
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

install: $(SIDL_FILES) $(MAKE_INC) install-recursive
	for file in $(SIDL_FILES); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(includedir)/$$file ; \
	done
	$(INSTALL_DATA) $(MAKE_INC) $(DESTDIR)$(libdir)/$(MAKE_INC)
	$(ECHO) "iGeom_SIDL_includes=-I$(includedir) -I$(includedir)/SIDL" >> $(DESTDIR)$(libdir)/$(MAKE_INC)
	$(ECHO) "iGeom_SIDL_LFDFLAGS=-L${libdir}" >> $(DESTDIR)$(libdir)/$(MAKE_INC)

uninstall: uninstall-recursive
	for file in $(SIDL_FILES); do \
	  rm -f $(DESTDIR)$(includedir)/$$file ; \
	done
	rm -f $(DESTDIR)$(libdir)/$(MAKE_INC)

install-recursive: recursive
	if test "x$(USE_BABEL)" = "xyes"; then \
	  for lib in gserver Cclient Fclient; do \
	    (cd $$lib && \
	      $(MAKE) -f ../$(srcdir)/$${lib}.make install \
	         DESTDIR="$(DESTDIR)" \
	         prefix="$(prefix)" \
	         exec_prefix="$(exec_prefix)" \
	         libdir="$(libdir)" \
	         datarootdir="$(datarootdir)" \
	         datadir="$(datadir)" \
	         includedir="$(includedir)" \
	    ) || exit 1; \
	  done ; \
	fi

uninstall-recursive: recursive
	if test "x$(USE_BABEL)" = "xyes"; then \
	  for lib in gserver Cclient Fclient; do \
	    (cd $$lib && \
	      $(MAKE) -f ../$(srcdir)/$${lib}.make uninstall \
	         DESTDIR="$(DESTDIR)" \
	         prefix="$(prefix)" \
	         exec_prefix="$(exec_prefix)" \
	         libdir="$(libdir)" \
	         datarootdir="$(datarootdir)" \
	         datadir="$(datadir)" \
	         includedir="$(includedir)" \
	    ) || exit 1; \
	  done ; \
	fi

distdir: $(DIST_FILES)
	test -d $(distdir) || $(mkdir_p) $(distdir) || exit 1
	test -d $(distdir)/gserver || $(mkdir_p) $(distdir)/gserver || exit 1
	@for file in $(DIST_FILES); do \
	  $(ECHO) "cp -p $(srcdir)/$$file $(distdir)/$$file"; \
	  cp -p $(srcdir)/$$file $(distdir)/$$file || exit 1; \
	done

recursive: gserver/babel.make Cclient/babel.make Fclient/babel.make common.make

# The --vpath option is broken for some versions of babel, so
# copy files when doing out-of-tree build.
# These should do nothing for in-tree builds because the target 
# and source files are the same file and therefore always have
# the same timestamp.
# Use 'cat' rather than 'cp' so that target files get a clean
# set of permissions (writable, owned by person running make, etc.)
gserver/babel.make: $(SIDL_FILES_PATH) $(SERVER_FILES_PATH)
	test -d gserver || $(mkdir_p) gserver || exit 1
	@if test "x$(srcdir)" != "x."; then \
	  for p in $(SERVER_FILES); do \
	    $(ECHO) "cat $(srcdir)/$$p > $$p"; \
	    cat $(srcdir)/$$p > $$p || exit 1; \
	  done \
	fi
	@if test "x$(USE_BABEL)" = "xyes"; then \
	  $(ECHO) '$(BABEL) $(BABELFLAGS) -sC++ $(SIDL_FILES_PATH) -o gserver' ; \
	  $(BABEL) $(BABELFLAGS) -sC++ $(SIDL_FILES_PATH) -o gserver && touch $@; \
	fi

Cclient/babel.make: $(SIDL_FILES)
	test -d Cclient || $(mkdir_p) Cclient || exit 1
	@if test "x$(USE_BABEL)" = "xyes"; then \
	  $(ECHO) '$(BABEL) $(BABELFLAGS) -cC $(SIDL_FILES_PATH) -o Cclient' ; \
	  $(BABEL) $(BABELFLAGS) -cC $(SIDL_FILES_PATH) -o Cclient && touch $@; \
	fi

Fclient/babel.make: $(SIDL_FILES)
	test -d Fclient || $(mkdir_p) Fclient || exit 1
	@if test "x$(USE_BABEL)" = "xyes"; then \
	  $(ECHO) '$(BABEL) $(BABELFLAGS) -cF77 $(SIDL_FILES_PATH) -o Fclient' ; \
	  $(BABEL) $(BABELFLAGS) -cF77 $(SIDL_FILES_PATH) -o Fclient && touch $@; \
	fi

.SUFFIXES: .cpp .cc .c .lo

check: 
	@test "x$(USE_BABEL)" != "xyes" || $(MAKE) run_testgeom

run_testgeom: testgeom
	@$(ECHO) "./$<"
	@./$< || $(ECHO) "<<<<<<<<<<<<<<<< testgeom FAILED >>>>>>>>>>>>>>>>>>>"

testgeom : testgeom.lo gserver/libiGeomserver.la
	$(CXXLINK) $< gserver/libiGeomserver.la

gserver/libiGeomserver.la : all

testgeom.lo : gserver/babel.make

.cpp.lo:
	$(LTCXXCOMPILE) -c -o $@ $<

.c.lo:
	$(LTCOMPILE) -c -o $@ $<

dvi:

pdf:

ps:

html:

info:

man:

installcheck:

.PHONEY: all clean distclean mostlyclean distdir check recursive \
         install install-recursive uninstall uninstall-recursive \
         dvi ps pdf html info man installcheck



