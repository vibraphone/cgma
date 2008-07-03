vpath = gserver bserver 
VPATH = gserver bserver 
BABEL_FLAGS = -V$(srcdir)

LD_FLAGS_PIC = -fpic -shared
SERVER_HDRS = iBase_ArrTag_IOR.h iBase_CreationStatus_IOR.h iBase_EntSet_IOR.h    \
  iBase_EntTag_IOR.h iBase_EntityType_IOR.h iBase_ErrorActions_IOR.h          \
  iBase_ErrorType_IOR.h iBase_IOR.h iBase_SetBoolOps_IOR.h  \
  iBase_SetRelation_IOR.h iBase_SetTag_IOR.h iBase_StorageOrder_IOR.h         \
  iBase_TagValueType_IOR.h iBase_Tag_IOR.h iMesh_AdjacencyInfo_IOR.h          \
  iGeom_Construct_IOR.h iGeom_CoreQuery_IOR.h iGeom_Factory_IOR.h             \
  iGeom_Geometry_IOR.h iGeom_IOR.h iGeom_Iterators_IOR.h iGeom_Modify_IOR.h   \
  iGeom_Parametric_IOR.h iGeom_Primitives_IOR.h iGeom_SIDL_GeomSidl_IOR.h     \
  iGeom_SIDL_IOR.h iGeom_Shape_IOR.h iGeom_Tolerance_IOR.h                    \
  iGeom_Topology_IOR.h iGeom_Transforms_IOR.h \
  iBase.hh iBase_ArrTag.hh iBase_CreationStatus.hh iBase_EntSet.hh   \
  iBase_EntTag.hh iBase_EntityType.hh iBase_ErrorActions.hh    \
  iBase_ErrorType.hh iBase_SetBoolOps.hh iBase_SetRelation.hh iBase_SetTag.hh \
  iBase_StorageOrder.hh iBase_Tag.hh iBase_TagValueType.hh iMesh.hh           \
  iGeom_Booleans.hh iGeom_Construct.hh iGeom_CoreQuery.hh iGeom_Factory.hh    \
  iGeom_Geometry.hh iGeom_Iterators.hh iGeom_Modify.hh iGeom_Parametric.hh    \
  iGeom_Primitives.hh iGeom_SIDL.hh iGeom_SIDL_GeomSidl.hh iGeom_Shape.hh     \
  iGeom_Tolerance.hh iGeom_Topology.hh iGeom_Transforms.hh sidl.hh            \
  sidl_BaseClass.hh sidl_BaseException.hh sidl_BaseInterface.hh               \
  sidl_ClassInfo.hh sidl_ClassInfoI.hh sidl_DFinder.hh  \
  sidl_DLL.hh sidl_Finder.hh sidl_InvViolation.hh sidl_Loader.hh              \
  sidl_PostViolation.hh sidl_PreViolation.hh sidl_Resolve.hh                  \
  sidl_SIDLException.hh sidl_Scope.hh sidl_io.hh sidl_io_Deserializer.hh      \
  sidl_io_IOException.hh sidl_io_Serializeable.hh sidl_io_Serializer.hh       \
  sidl_rmi.hh sidl_rmi_ConnectRegistry.hh sidl_rmi_InstanceHandle.hh          \
  sidl_rmi_InstanceRegistry.hh sidl_rmi_Invocation.hh                         \
  sidl_rmi_NetworkException.hh sidl_rmi_ProtocolFactory.hh                    \
  sidl_rmi_Response.hh iGeom_SIDL_GeomSidl_Impl.hh iBase_Error_Impl.hh        \
  iBase_Error_IOR.h iBase_Error.hh 
SERVER_CXX_SRCS =                                   \
  iBase_ArrTag.cc iBase_EntSet.cc iBase_EntTag.cc                             \
  iBase_SetBoolOps.cc iBase_SetRelation.cc iBase_SetTag.cc iBase_Tag.cc       \
  iGeom_Booleans.cc iGeom_Construct.cc iGeom_CoreQuery.cc iGeom_Factory.cc    \
  iGeom_Geometry.cc iGeom_Iterators.cc iGeom_Modify.cc iGeom_Parametric.cc    \
  iGeom_Primitives.cc iGeom_SIDL_GeomSidl.cc iGeom_Shape.cc                   \
  iGeom_Tolerance.cc iGeom_Topology.cc iGeom_Transforms.cc sidl_BaseClass.cc  \
  sidle_BaseException.cc \
  sidl_BaseInterface.cc sidl_ClassInfo.cc sidl_ClassInfoI.cc sidl_DFinder.cc  \
  sidl_DLL.cc sidl_Finder.cc sidl_InvViolation.cc sidl_Loader.cc              \
  sidl_PostViolation.cc sidl_PreViolation.cc sidl_SIDLException.cc            \
  sidl_io_Deserializer.cc sidl_io_IOException.cc sidl_io_Serializeable.cc     \
  sidl_io_Serializer.cc sidl_rmi_ConnectRegistry.cc                           \
  sidl_rmi_InstanceHandle.cc sidl_rmi_InstanceRegistry.cc                     \
  sidl_rmi_Invocation.cc sidl_rmi_NetworkException.cc                         \
  sidl_rmi_ProtocolFactory.cc sidl_rmi_Response.cc iGeom_SIDL_GeomSidl_Skel.cc \
  iGeom_SIDL_GeomSidl_Impl.cc iBase_Error_Impl.cc \
  iBase_Error_Skel.cc iBase_Error.cc 
SERVER_C_SRCS = iGeom_SIDL_GeomSidl_IOR.c  iBase_Error_IOR.c 
	
SERVER_OFILES = $(SERVER_CXX_SRCES:.cc=.o) $(SERVER_C_SRCS:.c=.o)

#bin_PROGRAMS = testgeom
#testgeom_SOURCES = testgeom.cpp

LDADD = @top_srcdir@/itaps/libiGeom.la @abs_srcdir@/libiGeomserver.so ${BABEL_DIR}/lib/libsidl.la
#testgeom_DEPENDENCIES = @top_srcdir@/itaps/libiGeom.la @abs_srcdir@/libiGeomserver.so ${BABEL_DIR}/lib/libsidl.la

all: libiGeomserver.so

settings:
	@echo "FACTORY_OFILES = $(FACTORY_OFILES)"
	@echo "SERVER_OFILES  = $(SERVER_OFILES)"

#============= iBase =============
repo/.btimestamp: $(srcdir)/iBase.sidl
	-rm -rf repo
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -tXML -o repo $(srcdir)/iBase.sidl 
	touch $@ 

bserver/babel.make: repo/.btimestamp $(srcdir)/iBase.sidl
	-rm -f $@
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -sC++ -o bserver $(srcdir)/iBase.sidl

#============= iBase =============

repo/.gtimestamp: $(srcdir)/iGeom.sidl bserver/babel.make
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -tXML -o repo -Rrepo $(srcdir)/iGeom.sidl
	touch $@ 

gserver/babel.make: repo/.gtimestamp repo/.btimestamp $(srcdir)/iGeom_SIDL.sidl $(srcdir)/iGeom.sidl 
	-rm -f $@
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -R"repo" -sC++ -o gserver $(srcdir)/iGeom_SIDL.sidl

Cclient: repo/.gtimestamp repo/.gtimestamp $(srcdir)/iGeom_SIDL.sidl
	-rm -rf Cclient
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -R"repo" -cC -o Cclient $(srcdir)/iGeom_SIDL.sidl

F77client: repo/.gtimestamp repo/.btimestamp $(srcdir)/iGeom_SIDL.sidl
	-rm -rf F77client
	@BABEL_DIR@/bin/babel $(BABEL_FLAGS) -R"repo" -cF77 -o F77client $(srcdir)/iGeom_SIDL.sidl



libiGeomserver.so: gserver/babel.make ${SERVER_OFILES}
	$(CXX) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${SERVER_OFILES}

libiGeomCclient.so: Cclient ${CCLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${CCLIENT_OFILES}

libiGeomF77client.so: F77client ${F77CLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${F77CLIENT_OFILES}

Cclient/babel.make: Cclient

F77client/babel.make: F77client

clean:
	-rm -f *.o *.so testgeom *~ core.*

clean_all: clean
	-rm -rf F77client Cclient repo
	cd mserver; rm -f ${SERVER_HDRS} ${SERVER_SRCS}

.cxx.o:
	$(CXX) $(LD_FLAGS_PIC) $(CXX_FLAGS) -I$(BABEL_DIR)/include ${iGeom_SIDL_INCLUDES} -o $@ -c $<

.cpp.o:
	$(CXX) $(LD_FLAGS_PIC) $(CXX_FLAGS) -I$(BABEL_DIR)/include ${iGeom_SIDL_INCLUDES} -o $@ -c $<

.cc.o:
	$(CXX) $(LD_FLAGS_PIC) $(CXX_FLAGS) -I$(BABEL_DIR)/include ${iGeom_SIDL_INCLUDES} -o $@ -c $<

.c.o:
	$(CC) $(LD_FLAGS_PIC) $(CC_FLAGS) -I$(BABEL_DIR)/include ${iGeom_SIDL_INCLUDES} -o $@ -c $<

.f.o:
	$(FC) $(LD_FLAGS_PIC) $(FC_FLAGS) ${iGeom_SIDL_INCLUDES} -o $@ -c $<

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iGeom-SIDL-Defs.inc

# By default, iGeom-SIDL-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
iGeom-SIDL-Defs.inc:
	@echo "include $(cfgdir)/iGeom-SIDL-Defs.inc" > $@
	@echo 'iGeom_SIDL_INCLUDES = $$(iGeom_INCLUDES) -I$(BABEL_DIR)/include' >> $@
	@echo 'iGeom_SIDL_CPPFLAGS = $$(iGeom_CPPFLAGS) -I$(BABEL_DIR)/include' >> $@
	@echo 'iGeom_SIDL_DEFINES  = $$(iGeom_DEFINES)' >> $@
	@echo 'iGeom_SIDL_LDFLAGS  = $$(iGeom_LDFLAGS) -L$(BABEL_DIR)/lib' >> $@
	@echo 'iGeom_SIDL_LTFLAGS  = $$(iGeom_LDFLAGS) -R$(BABEL_DIR)/lib' >> $@
	@echo 'iGeom_SIDL_LIBS     = $$(iGeom_LIBS) -lsidl -liGeomserver' >> $@
