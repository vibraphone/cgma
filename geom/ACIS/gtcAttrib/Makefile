##---------------------------------------------------------------------------
## $Id: Makefile,v 1.1 2006/06/09 13:38:52 rakerr Exp $
##===========================================================================

SOURCE_DIR = ./src

LIB_NAME = gtcAttrib
LEAF_INCL = ./incl
LIB_INCL = -I$(LEAF_INCL)

SOURCE_FILES = \
attrib_edge_parametric_length.C \
attrib_gtc.C \
attrib_mark.C \
attrib_name.C 

SOURCE = $(SOURCE_FILES:%=$(SOURCE_DIR)/%)

VERSION = 1.0

ACIS_TOP_DIR = /usr/local/acis
AARCH = $(MACHINE_TYPE)

GTO_VERSION = 3.9
GTODIR = /local/workspace/goodyear/gto_$(GTO_VERSION)

GTOINCLUDEDIR = $(GTODIR)/include/$(GTO_VERSION)

GTOLIBDIR = $(GTODIR)/lib/$(MACHINE_TYPE)

include .tce.makefile.acis.R16
MAKEDEPEND = makedepend
include Makefile.$(MACHINE_TYPE)

OBJECT_DIR = ./obj_$(MACHINE_TYPE)
SHARED_OBJECT_DIR = $(OBJECT_DIR)/shrd


DEVELOPER = 

INCL     = $(LIB_INCL) -I$(STL_HOME) $(ACIS_INCLUDES)


OBJ = $(SOURCE:$(SOURCE_DIR)/%.C=$(OBJECT_DIR)/%.o)
SHARED_OBJ = $(SOURCE:$(SOURCE_DIR)/%.C=$(SHARED_OBJECT_DIR)/%.o)

.SUFFIXES: .C .h .H .o .i

objs : $(OBJ)
shrd_objs : $(SHARED_OBJ)

#this is a cheesy way to see what'ts being macro expanded.
#.C.o:
#	$(CCC) -c -E $(CXXFLAG) $(COMPILEC) $(INCL) $< > $@.E

%.o:$(SOURCE_DIR)/%.C
	${MAKE} $(OBJECT_DIR)/$*.o

$(OBJECT_DIR)/%.o:$(SOURCE_DIR)/%.C
	@echo ' '
	@echo Making $@
	@echo ' '
	$(CCC) -c $(CXXFLAG) $(UCXXFLAG) $(DEVELOPER) $(COMPILEC) $(INCL) -o $@ $<

$(SHARED_OBJECT_DIR)/%.o:$(SOURCE_DIR)/%.C
	$(MAKE) UCXXFLAG='$(SHARED_UCFLAG)' OBJECT_DIR='$(SHARED_OBJECT_DIR)' $@

$(LIB_NAME)_exConst.o:$(OBJECT_DIR)/$(LIB_NAME)_exConst.o

$(OBJECT_DIR)/$(LIB_NAME)_exConst.o:$(SOURCE_DIR)/$(LIB_NAME)_exConst.C
	$(CCC) -c $(CXXFLAG) $(UCXXFLAG) $(DEVELOPER) $(COMPILEC) $(INCL) -o $@ $<

$(SOURCE_DIR)/$(LIB_NAME)_exConst.C : $(LEAF_INCL)/$(LIB_NAME)_exConst.h
	rm -f $@
	$(ECHO) '#ifndef $(LIB_NAME)_EXCONST_C__\n#define $(LIB_NAME)_EXCONST_C__\n\n#include <$(LIB_NAME)_exConst.h>\n\nstring $(LIB_NAME)_version = "'$(VERSION)'" ;\nstring $(LIB_NAME)_acis_v = "'$(ACIS_VERSION)'" ;\nstring $(LIB_NAME)_date = "'`date`'" ;\n\nstring $(LIB_NAME)_gto_version = "'$(GTO_VERSION)'" ;\nstring $(LIB_NAME)_compiler_string = "$(COMPILER_STRING)" ;\n\n#endif\n'  > $@


%.i:$(SOURCE_DIR)/%.C 
	$(CCC) -c -E $(CXXFLAG) $(UCXXFLAG) $(DEVELOPER) $(COMPILEC) $(INCL) $< -o $@

lib: lib/$(MACHINE_TYPE)/lib$(LIB_NAME)$(VERSION).a

lib_so: lib/$(MACHINE_TYPE)/lib$(LIB_NAME)$(VERSION).so

lib/$(MACHINE_TYPE):
	mkdir -p lib/$(MACHINE_TYPE)

lib/$(MACHINE_TYPE)/lib$(LIB_NAME)$(VERSION).a : $(OBJ) lib/$(MACHINE_TYPE)
	touch $(LEAF_INCL)/$(LIB_NAME)_exConst.h
	$(MAKE) $(LIB_NAME)_exConst.o
	rm -f $@
	$(ARCHIVER) $(LIBFLAGS) $@ $(OBJ) $(OBJECT_DIR)/$(LIB_NAME)_exConst.o $(ARCH_OBJ)

lib/$(MACHINE_TYPE)/lib$(LIB_NAME)$(VERSION).so : lib/$(MACHINE_TYPE) FORCE
	$(MAKE) clean
	$(MAKE) UCXXFLAG='$(SHARED_UCFLAG)' OBJECT_DIR='$(SHARED_OBJECT_DIR)'
	touch $(LEAF_INCL)/$(LIB_NAME)_exConst.h
	$(MAKE) UCXXFLAG='$(SHARED_UCFLAG)' OBJECT_DIR='$(SHARED_OBJECT_DIR)' $(LIB_NAME)_exConst.o
	rm -f $@
	$(SHARED_ARCHIVER) $(SHARED_LFLAGS) $@ $(SHARED_OBJ) $(SHARED_OBJECT_DIR)/$(LIB_NAME)_exConst.o

install::
	$(MAKE) UCXXFLAG='$(LIB_UCXXFLAG)' DEVELOPER='-DNDEBUG' lib
	$(MAKE) DEVELOPER='-DNDEBUG' lib_so

clean: 
	rm -rf $(OBJ) $(SHARED_OBJ) $(ARCH_CLEAN) 

cleansource: 
	rm -rf incl/*~ src/*~

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

depend: make.depend.$(MACHINE_TYPE)

make.depend.$(MACHINE_TYPE):$(OBJECT_DIR)
	touch $@
	${MAKEDEPEND} -f $@ $(SOURCE) -Y -I$(LEAF_INCL) -I$(GTOINCLUDEDIR) -I$(GTO2ACIS_DIR)/incl $(LIB_INCL) 1>>/dev/null 2>&1
	sed 's^$(SOURCE_DIR)^$(OBJECT_DIR)^g' $@ > $@.tmp
	mv $@.tmp $@

distrib:
	$(MAKE) clean

FORCE:



include make.depend.$(MACHINE_TYPE)

##---------------------------------------------------------------------------
## $Log: Makefile,v $
## Revision 1.1  2006/06/09 13:38:52  rakerr
## Adding gtcAttrib files.
##
## Revision 1.2  2006/05/30 14:16:24  prwolfe
## New version numbners and added check for get_name not to return the generic
##    names managed by the gto backpointer attribute.
##
## Revision 1.1  2006/05/03 15:15:24  prwolfe
## Initial setup and merge of file maintained by gto2acis, cubit, and qmesh.
##    this library will give us a single point of maintenence and argument.
##
## Revision 1.4  2001/11/28 21:43:54  prwolfe
## Changed from ARCH to MACHINE_TYPE - prw.
##
## Revision 1.3  2001/01/19 15:33:05  a314773
## New source file full_depth_sipes.C - prw.
##
## Revision 1.2  2000/08/25 11:23:57  a314773
## Added blind_hole_sipe_surf.C to the sources - dkk.
##
## Revision 1.1  2000/03/23 10:59:43  a314773
## Makefile for the threeD library.  Sets a few items and then includes
## 	Makefile.leaf from the top level directory.
##
##
##===========================================================================
