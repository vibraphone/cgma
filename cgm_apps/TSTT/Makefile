VPATH = server Cclient F77client
vpath = server Cclient F77client

include TSTTG-Defs.inc


######## include make.rules
%.d : %.cxx
	$(MKDEP) -DIS_BUILDING_MB $(MKDEP_FLAGS) $(CXXFLAGS) $< > /dev/null

%.d : %.cpp
	$(MKDEP) -DIS_BUILDING_MB $(MKDEP_FLAGS) $(CXXFLAGS) $< > /dev/null

%.d : %.cc
	$(MKDEP) -DIS_BUILDING_MB $(MKDEP_FLAGS) $(CXX_FLAGS) $< > /dev/null

%.d : %.c
	$(MKDEP) -DIS_BUILDING_MB $(MKDEP_FLAGS) $(CC_FLAGS) $< > /dev/null

%.d : %.f
	$(MKDEP) -DIS_BUILDING_MB $(MKDEP_FLAGS) $(FC_FLAGS) $< > /dev/null

%.o : %.cxx
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.cpp
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.cc
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.c
	$(CC) -DIS_BUILDING_MB $(CC_FLAGS) -o $@ -c $<

%.o : %.f
	$(FC) -DIS_BUILDING_MB $(FC_FLAGS) -o $@ -c $<

###############################################################

######## include make.Linux
CXX = g++
CXX_FLAGS_REQ = -fpic
CXX_FLAGS_DBG = -g -Wall
CXX_FLAGS_OPT = -O3 -funroll-loops -fexpensive-optimizations

CC  = gcc
CC_FLAGS_REQ = -fpic
CC_FLAGS_DBG = -g -Wall
CC_FLAGS_OPT = -O3

FC  = g77
#FC   = /usr/apps/pgi/3.3/linux86/bin/pgf77
FC_FLAGS_REQ = 
FC_FLAGS_DBG = -g
FC_FLAGS_OPT = -O

MKDEP = gcc -E -MMD

LD = g++
LD_FLAGS_REQ = -fpic -shared
LD_FLAGS_DBG = -g
LD_FLAGS_OPT =

ARCHIVER = ar cr


INCPATH = ${TSTTG_SERVER_INCLUDES}

CXX_FLAGS = $(CXX_FLAGS_REQ) $(CXX_FLAGS_DBG) $(INCPATH) $(MACH_CXXFLAGS)
CC_FLAGS = $(CC_FLAGS_REQ) $(CC_FLAGS_DBG) $(INCPATH)
FC_FLAGS = $(FC_FLAGS_REQ) $(FC_FLAGS_DBG) $(INCPATH)
LD_FLAGS = $(LD_FLAGS_REQ) $(LD_FLAGS_DBG) $(INCPATH)
###########################################################

include server/babel.make
SERVER_IMPL_HDRS := $(subst TSTTB_Error_Impl.hh,,$(IMPLHDRS))
SERVER_IMPL_SRCS := $(subst TSTTB_Error_Impl.cc,,$(IMPLSRCS))
SERVER_GEN_IOR_HDRS := $(subst TSTTB_Error_IOR.h,,$(IORHDRS))
SERVER_GEN_IOR_SRCS := $(subst TSTTB_Error_IOR.c,,$(IORSRCS))
SERVER_GEN_STUB_HDRS := $(subst TSTTB_Error.hh,,$(STUBHDRS))
SERVER_GEN_STUB_SRCS := $(subst TSTTB_Error.cc,,$(STUBSRCS))
SERVER_GEN_SKEL_SRCS := $(subst TSTTB_Error_Skel.cc,,$(SKELSRCS))
SERVER_GEN_HDRS := $(SERVER_GEN_IOR_HDRS) $(SERVER_GEN_STUB_HDRS)
SERVER_GEN_SRCS := $(SERVER_GEN_IOR_SRCS) $(SERVER_GEN_STUB_SRCS) $(SERVER_GEN_SKEL_SRCS) 
SERVER_OFILES := $(subst .cc,.o,$(SERVER_GEN_SRCS) $(SERVER_IMPL_SRCS))
SERVER_OFILES := $(subst .c,.o, $(SERVER_OFILES))

-include Cclient/babel.make
CCLIENT_GEN_HDRS = $(IORHDRS) $(STUBHDRS)
CCLIENT_GEN_SRCS = $(STUBSRCS)
CCLIENT_OFILES := $(subst .cc,.o,$(CCLIENT_GEN_SRCS))
CCLIENT_OFILES := $(subst .c,.o, $(CCLIENT_OFILES))

-include F77client/babel.make
F77CLIENT_GEN_HDRS = $(IORHDRS) $(STUBHDRS)
F77CLIENT_GEN_SRCS = $(STUBSRCS)
F77CLIENT_OFILES := $(subst .cc,.o,$(F77CLIENT_GEN_SRCS))
F77CLIENT_OFILES := $(subst .c,.o, $(F77CLIENT_OFILES))

all: testgeom

repo/.timestamp: TSTTG.sidl
	$(BABEL) -tXML -R${TSTTB_DIR}/repo -o repo TSTTG.sidl
	touch $@ 

server/babel.make: repo/.timestamp ${TSTTB_DIR}/repo/.timestamp TSTTG_CGM.sidl TSTTG.sidl
	-rm -f $@
	$(BABEL) -R"repo;${TSTTB_DIR}/repo" -sC++ -o server TSTTG_CGM.sidl

Cclient: repo/.timestamp ${TSTTB_DIR}/repo/.timestamp TSTTG_CGM.sidl
	-rm -rf Cclient
	$(BABEL) -R"repo;${TSTTB_DIR}/repo" -cC -o Cclient TSTTG_CGM.sidl

F77client: repo/.timestamp ${TSTTB_DIR}/repo/.timestamp TSTTG_CGM.sidl
	-rm -rf F77client
	$(BABEL) -R"repo;${TSTTB_DIR}/repo" -cF77 -o F77client TSTTG_CGM.sidl

libCGMserver.so: server/babel.make ${SERVER_OFILES} TSTTG_CGM.o CATSTT.o server/babel.make
	$(LD) $(LD_FLAGS) -o $@ ${SERVER_OFILES} TSTTG_CGM.o CATSTT.o

libCGMserver.a: server/babel.make ${SERVER_OFILES} TSTTG_CGM.o CATSTT.o server/babel.make
	$(ARCHIVER) $@ ${SERVER_OFILES} TSTTG_CGM.o CATSTT.o

libCGMCclient.so: Cclient ${CCLIENT_OFILES}
	$(LD) $(LD_FLAGS) -o $@ ${CCLIENT_OFILES}

libCGMF77client.so: F77client ${F77CLIENT_OFILES}
	$(LD) $(LD_FLAGS) -o $@ ${F77CLIENT_OFILES}

Cclient/babel.make: Cclient

F77client/babel.make: F77client


testgeom: testgeom.o ${TSTTG_SERVER_FILES}
	$(LD)  -g -o $@ testgeom.o -Wl,-rpath,. ${TSTTG_SERVER_LIBS}

testc: testc.o ${TSTTG_CLIENT_C_FILES} ${TSTTG_SERVER_FILES} ${TSTTB_SERVER_FILES} 
	$(LD)  -o $@ testc.o -Wl,-rpath,. ${TSTTG_CLIENT_C_LIBS} ${TSTTG_SERVER_LIBS}

clean:
	-rm -f *.o *.so testc testgeom *~ core.*

clean_all: clean
	-rm -rf F77client Cclient repo
	cd server; rm -f ${SERVER_GEN_HDRS} ${SERVER_GEN_SRCS}

depend:
	${MAKEDEPEND} ${INCPATH} testc.c testgeom.cpp TSTTG_CGM.cpp CATSTT.cpp > make.dependencies


-include make.dependencies
