#
# This file contains make macro definitions needed for a cubit developers
# installation.  This file is included directly into the CUBIT makefile,
# so definitions are in make macro format
#
CUBITROOT = $(CUBIT)
#
# SOURCE_DIR:
# Developer's primary cubit source directory (i.e. the directory from which
# cvs operations are done)
#
SOURCE_DIR = .

#
# X11_LIB_DIR, X11_INCLUDE:
# Where to look for X11 libs and includes 
#

X11_LIB_DIR = /usr/openwin/lib
X11_INCLUDE = -I/usr/openwin/include
#####For dynamic linking of libraries do
X11_LINK =  -L${X11_LIB_DIR} -lXt -lXext -lX11
#####For static linking do
#X11_LINK = -L${X11_LIB_DIR} -lXt -lX11 -Bstatic -lXext -Bdynamic
#
# EXODUS_LIB_DIR, EXODUS_INCLUDE:
# Where to look for EXODUS includes and libs
#
EXODUS_INCLUDE =   -I$(CUBITROOT)/exodus/exodus-4.01/include
EXODUS_LIB_DIR = $(CUBITROOT)/exodus/exodus-4.01/lib
EXODUS_LIBS = -lexoIIv2c401

NETCDF_INCLUDE = -I$(CUBITROOT)/netcdf/netcdf-3.4.snl3/include
NETCDF_LIB_DIR = $(CUBITROOT)/netcdf/netcdf-3.4.snl3/lib
NETCDF_LIBS = -lnetcdf

# Graph Drawing is used for viewing the Whisker weaving process
# and is not necessary for all of CUBIT.  If you do not want it,
# Then you can take the link out of the makefile (GRAPH_DRAWING_LINK).
# Also you will need to take it out of the defines section in
# the MACH_CXXFLAGS
#GRAPH_DRAWING_INCLUDE = -I$(CUBITROOT)/include 
#GRAPH_DRAWING_LIB_DIR = $(CUBITROOT)/lib 
#GRAPH_DRAWING_LINK = -L${GRAPH_DRAWING_LIB_DIR} -lgraphpack

#
# LP_SOLVE_INCLUDE and LP_SOLVE_LIB_DIR
# Where to look for LP Solve includes and libs
#
LP_SOLVE_INCLUDE = -I$(CUBITROOT)/include
LP_SOLVE_LIB_DIR = $(CUBITROOT)/lib

##Workshop 7
#SUN_WORKSHOP=/net/sasn232-atm/opt/SUNWspro7.0/SUNWspro
#CC_INCLUDE = -I$(SUN_WORKSHOP)/prod/include/CCios \
#	-I$(SUN_WORKSHOP)/prod/include/CC/Cstd \
#	-I$(SUN_WORKSHOP)/prod/include/CC

##Studio 8
SUN_WORKSHOP=/opt/Studio8/SUNWspro
CC_INCLUDE = -I$(SUN_WORKSHOP)/prod/include/CCios \
	-I$(SUN_WORKSHOP)/prod/include/CC/Cstd \
	-I$(SUN_WORKSHOP)/prod/include/CC


SYS_INCLUDE = -I/usr/include/sys

# Sometimes if we put an include path in the previous variable,
# the compile goes haywire.  But if we don't put the paths somewhere
# then makedepend can't find them.  This variable is kind of a dummy
# so that we can tell makedepend where to find things, but not have
# it affect the actual compile.
MAKE_DEPEND_INCLUDE = 
 
# These include macros are for basic cubit versions. To set for
# parallel operations (MPI), create a COMM_INCLUDE_OPTION
# with the definition, -I$(MPI_INCLUDE_DIR)
 
EXTRA_INCLUDE = $(COMM_INCLUDE_OPTION) -I/usr/include

#
# OptMS Package
#
# OPTMS_OPTION should be ON or OFF, indicating whether to link
# in the library or not.
#
#   **IMPORTANT!!!**
#  After changing the value of OPTMS_OPTION, you need to type
#  'make optms_files', or the affected files won't get recompiled!!!!
#
OPTMS_OPTION = OFF
OPTMS_ROOT = $(CUBITROOT)/Opt-MS
OPTMS_INCLUDE = -I$(OPTMS_ROOT)/include
#FORTRAN_LIB =  -L$(SUN_WORKSHOP)/SC4.2/lib -lF77 -lM77 -lsunmath 

#MESQUITE OPTIONS
MESQUITE_VER = 08c1
MESQUITE_DIR = $(CUBITROOT)/mesquite/mesquite_$(MESQUITE_VER)
MESQUITE_OPTION = ON
MESQUITE_INCLUDE = -I$(MESQUITE_DIR)/include
MESQUITE_LIB_DIR = $(MESQUITE_DIR)/lib/Sun
MESQUITE_LIB     = -lmesquite$(MESQUITE_VER)
MESQUITE_FLAGS = -DUSE_C_PREFIX_INCLUDES -DUSE_STD_INCLUDES

#Workshop 7 (Fortran 90)
#static
#FORTRAN_LIB =	-Bstatic -lfsu -lfui -lf77compat -lsunmath -Bdynamic
#dynamic
FORTRAN_LIB =	-lfsu -lfui -lf77compat -lmvec -lsunmath
APP_LIB = 

OPTMS_LINK = $(OPTMS_ROOT)/lib/libO/solaris/libOPTMS.a \
             $(OPTMS_ROOT)/lib/libO/solaris/blas.a \
            $(FORTRAN_LIB) $(APP_LIB)

# camal variables
CAMAL_LIB_DIR	= $(SOURCE_DIR)/camal/lib/SunOS

# Simulog/INRIA tet sources; set to ON to link in INRIA tet mesher,
# otherwise set to OFF
SIMULOG_OPTION = ON
SIMULOG_LIB_DIR = $(CUBITROOT)/camal/camal2.0.2/lib/SunOS

# feature-based decomp'n sources; set to ON to link in feature-based decomp'n,
# otherwise set to OFF
FEATURE_OPTION = ON

# Medial sources; set to ON to link in the medial axis based work.
# otherwise set to OFF
MEDIAL_OPTION = OFF
MEDIAL_DEFINE =
MEDIAL_INCLUDES =
#MEDIAL_OPTION = ON
#MEDIAL_DEFINE = -DUSING_MEDIAL
#MEDIAL_INCLUDES = -Imedial -Imedial/medial_util -Imedial/medial2D \
#	-Imedial/medial3D


# verdict variables
CUBIT_VERDICT_VERSION = 112
VERDICT_DIR = $(CUBITROOT)/verdict/verdict1.1.2
VERDICT_LIB_DIR =  $(VERDICT_DIR)/lib
VERDICT_LIB = -lverdict$(CUBIT_VERDICT_VERSION)
VERDICT_INCLUDE =  -I${VERDICT_DIR}/include


# showviz settings
VTK_DIR = $(CUBITROOT)/VTK/VTK-4.2.2
VTK_LINK = -L${VTK_DIR}/lib/vtk -lvtkHybrid -lvtkIO -lvtkRendering -lvtkftgl -lvtkfreetype\
           -lvtkImaging -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtkexpat \
           -lvtkjpeg -lvtkpng -lvtktiff -lvtkzlib -lGL -lpthread -lrt
VTK_INCLUDE = -I${VTK_DIR}/include/vtk -I$(SOURCE_DIR)/showviz/geom -I$(SOURCE_DIR)/showviz/mesh

SHOWVIZ_OPTION=ON
HOOPS_OPTION=OFF

# HOOPS variables - these variables are more system dependent than others -
# make sure to check these carefully 
#
# HOOPS needs additional system libraries linked in, which vary depending
# which system you're on
# for the SUN:
EXTRA_HOOPS_LIB =
CUBIT_HOOPS_VERSION = 620
HOOPS_DIR = $(CUBITROOT)/hoops/hoops-6.20
HOOPS_LIB_DIR =  $(HOOPS_DIR)/lib
HOOPS_INCLUDE =  -I$(HOOPS_DIR)/include
#For a static link use the next line
#HOOPS_LIB = -Bstatic -lhoops$(CUBIT_HOOPS_VERSION) -Bdynamic -ldl
#For a dynamic link use the next line
HOOPS_LIB = -lhoops610

EXTRA_LIBS_LINK =   

# ACIS directories
#
FACETER = -lfaceter

####Define which solid modeling engine you will use, or all of them...
SOLID_MODELER_DEFINES = -DACIS_3D -DACIS_HEALER -DACIS_LOCAL_OPS \
			-DACIS_IGES_TRANSLATOR \
			-DACIS_STEP_TRANSLATOR 

### Note that you can also use these, if you have the right license
# and are using acis 7.0.x
#-DACIS_CATIA_TRANSLATOR -DACIS_PROE_TRANSLATOR
# if you add them, then add the following to the ACIS_LINK below
#		-lcathusk \
#		-lproehusk \
#		-lxgeometric \
#

#
#ACIS 14.1 version
#
CUBIT_ACIS_VERSION = 1401
ACIS_VERSION = 1401
ACIS_DIR =  $(CUBITROOT)/acis/acis14.1
ACIS_SYSTEM = solaris_so
ACIS_LIB_DIR =  $(ACIS_DIR)/bin/$(ACIS_SYSTEM)

ACIS_INCLUDE = -I$(ACIS_DIR)/include

ACIS_ARCH = solaris

ACIS_LINK       = -L${ACIS_LIB_DIR}\
		-lacisstep\
		-lacisiges\
		-lxiges\
		-lxstep\
		-lxacis2k\
		-lxcore2k\
		-lSPAXAcisBase\
                -lSPAXAssemblyRep\
                -lSPAXDefaultGeometryRep\
		-lSpaAVis\
		-lSpaAWarp\
		-lSpaASurf\
		-lSpaALops\
		-lSpaABlend\
		-lSPAXInterop\
		-lSPAXBase\
		-lSpaACIS\
		-lSpaBase 


#
#ACIS 13.2 version
#
#CUBIT_ACIS_VERSION = 1302
#ACIS_DIR =  $(CUBITROOT)/acis/acis13.2
#ACIS_VERSION = 1302
#ACIS_SYSTEM = solaris_so
#ACIS_LIB_DIR =  $(ACIS_DIR)/bin/$(ACIS_SYSTEM)

#ACIS_INCLUDE = -I$(ACIS_DIR)/include

#ACIS_ARCH = solaris

#ACIS_LINK       = -L${ACIS_LIB_DIR}\
#		-lxstep\
#		-lxacis2k\
#		-lxcore2k\
#		-lxutil\
#		-lacisstep\
#		-lxiges\
#		-lacisiges\
#		-lSpaAWarp\
#		-lSpaAVis\
#		-lSpaASurf\
#		-lSpaALops\
#		-lSpaABlend\
#		-lSpaACIS\
#		-lSpaBase

#
#ACIS 7.0.5 version
#
#CUBIT ACIS version 7.0
#CUBIT_ACIS_VERSION = 700
#ACIS_DIR =  $(CUBITROOT)/acis/acis7.0.5
#ACIS_VERSION = 700
#ACIS_SYSTEM = solaris_so
#ACIS_LIB_DIR =  $(ACIS_DIR)/lib/$(ACIS_SYSTEM)

#ACIS_INCLUDE = -I$(ACIS_DIR) \
#		-I$(ACIS_DIR)/kern \
#		-I$(ACIS_DIR)/base \
#		-I$(ACIS_DIR)/law \
#		-I$(ACIS_DIR)/fct \
#		-I$(ACIS_DIR)/intr \
#		-I$(ACIS_DIR)/swp \
#		-I$(ACIS_DIR)/blnd \
#		-I$(ACIS_DIR)/covr \
#		-I$(ACIS_DIR)/bool \
#		-I$(ACIS_DIR)/cstr \
#		-I$(ACIS_DIR)/eulr \
#		-I$(ACIS_DIR)/ofst \
#		-I$(ACIS_DIR)/ga \
#		-I$(ACIS_DIR)/heal \
#		-I$(ACIS_DIR)/iges \
#		-I$(ACIS_DIR)/step \
#		-I$(ACIS_DIR)/lopt \
#		-I$(ACIS_DIR)/rem \
#		-I$(ACIS_DIR)/skin \
#		-I$(ACIS_DIR)/lop \
#		-I$(ACIS_DIR)/proe \
#		-I$(ACIS_DIR)/catia \
#		-I$(ACIS_DIR)/oper \
#		-I$(ACIS_DIR)/ct \
#		-I$(ACIS_DIR)/sbool 

#ACIS_ARCH = solaris
#ACIS_LINK       = -L${ACIS_LIB_DIR}\
#		-lhealhusk \
#		-ligeshusk \
#		-lstephusk \
#		-lcaselib \
#		-lshl_husk \
#		-llop_husk \
#		-lblend \
#		-lsweep \
#		-lsbool \
#		-lskin \
#		-lcover \
#		-loffset \
#		-ltransutl \
#		-loperator \
#		-lrem_husk \
#		-lrbi_husk \
#		-llopt_husk \
#		-lboolean \
#		-lskin \
#		-lga_husk \
#		-lrnd_husk \
#		-lct_husk \
#		-leuler \
#		-lconstrct \
#		-lfaceter \
#		-lintersct \
#		-lkernel \
#		-llawutil \
#		-lbaseutil 

#
# ACIS 6.x version
#
#SOLID_MODELER_DEFINES = -DACIS_3D -DACIS_HEALER -DACIS_LOCAL_OPS \
#		-DACIS_IGES_TRANSLATOR -DACIS_STEP_TRANSLATOR -DMMGR_FREELIST
#ACIS_VERSION = 600
#ACIS_SYSTEM = solaris_so

##Minor version dependent stuff
#### Acis 6.2 Stuff           ###
#CUBIT_ACIS_VERSION = 602
#ACIS_DIR =  $(CUBITROOT)/acis/acis6.2

#### Acis 6.1 Stuff           ###
#CUBIT_ACIS_VERSION = 601
#ACIS_DIR =  $(CUBITROOT)/acis/acis6.1

### Acis 6.0 Stuff           ###
#CUBIT_ACIS_VERSION = 600
#ACIS_DIR =  $(CUBITROOT)/acis/acis6.0


##Minor version INdependent stuff
##
######If you want shared libraries, do "solaris_so", else do "solaris"
### Shared libraries don't link in the libraries statically to the
### cubit executable.  This decreases executable size if you have
### acis running on the machines you will run CUBIT.

#ACIS_LIB_DIR =  $(ACIS_DIR)/lib/$(ACIS_SYSTEM)
##
#ACIS_INCLUDE = -I$(ACIS_DIR) \
#		-I$(ACIS_DIR)/kern \
#		-I$(ACIS_DIR)/base \
#		-I$(ACIS_DIR)/law \
#		-I$(ACIS_DIR)/fct \
#		-I$(ACIS_DIR)/intr \
#		-I$(ACIS_DIR)/swp \
#		-I$(ACIS_DIR)/blnd \
#		-I$(ACIS_DIR)/covr \
#		-I$(ACIS_DIR)/bool \
#		-I$(ACIS_DIR)/cstr \
#		-I$(ACIS_DIR)/eulr \
#                -I$(ACIS_DIR)/ofst \
#		-I$(ACIS_DIR)/ga \
#		-I$(ACIS_DIR)/heal \
#		-I$(ACIS_DIR)/iges \
#		-I$(ACIS_DIR)/step \
#		-I$(ACIS_DIR)/lopt \
#		-I$(ACIS_DIR)/rem \
#		-I$(ACIS_DIR)/skin \
#		-I$(ACIS_DIR)/lop \
#		-I$(ACIS_DIR)/oper \
#		-I$(ACIS_DIR)/mmgr

#ACIS_ARCH = solaris
#ACIS_LINK       = -L${ACIS_LIB_DIR}\
#		-lhealhusk \
#		-lblend \
#		-ligeshusk \
#		-lstephusk \
#		-lcaselib \
#		-llop_husk \
#		-lsweep \
#		-lsbool \
#		-lskin \
#		-lcover \
#		-loffset \
##		-ltransutl \
#		-loperator \
#		-lrem_husk \
#		-lrbi_husk \
#		-llopt_husk \
#		-lboolean \
#		-lskin \
#		-lga_husk \
#		-lrnd_husk \
#		-lct_husk \
#		-leuler \
#		-lconstrct \
#		-lfaceter \
#		-lintersct \
#		-lkernel \
#		-lspline \
#		-lkernel \
#		-llawutil \
#		-lbaseutil 


ACIS_DEFINES =  -D$(ACIS_ARCH) -D$(ACIS_SYSTEM) -DCXX30\
                -DACIS$(CUBIT_ACIS_VERSION) \
                -DCUBIT_ACIS_VERSION=$(CUBIT_ACIS_VERSION) -DPOINT_VERSION=0 \
                -DAGVER=105 -DANSI -DNO_DTORS

# If you want to optimise the code 
# DEBUG_FLAG = -fast -DNDEBUG -fsimple=1 -xtarget=generic
DEBUG_FLAG = -g

# Template files that need to be compiled differently depending on the 
# platform.  The Linux platform has to #include the .cpp file into the .hpp
# file, and so has to have it taken out of the makefile.  So, we declare the
# *_TEMPLATES variable and set it or not, depending on the platform.
LIST_TEMPLATES = DLIList.cpp \
		RTree.cpp \
		RTreeNode.cpp
FEATURE_TEMPLATES = FeatureTable.cpp
GEOM_TEMPLATES = BoundaryConstrainTool.cpp

# These compile flags are for basic cubit versions. To set for no graphics
# or parallel operations (MPI), add the -DNO_GRAPHICS and
# -DMPI_COMMUNICATIONS_ENABLED
# flags to GRAPHICS_OPTION and COMM_OPTION respectively

MACH_CFLAGS = $(ACIS_DEFINES) -DSOLARIS -DSUN -DSVR4 \
              -DSYSV -DANSI -DXTFUNCPROTO \
              $(GRAPHICS_OPTION) $(COMM_OPTION) -Dcplusplus \
              -DCUBIT_HOOPS_VERSION=$(CUBIT_HOOPS_VERSION)
 
# The -xpg is for profiling data, and if it isn't at the end of this
# and you want to run gprof you'll need to add it to the end of this 
# next line.
# If you want graph drawing uncomment, the first CXXFLAGS and comment
# the second.
#MACH_CXXFLAGS = $(ACIS_DEFINES) -D__EXTERN_C__ -DSOLARIS \
#                -DSUN -DSVR4 -DSYSV -DANSI -DXTFUNCPROTO +w \
#                -DGRAPH_DRAWING $(GRAPHICS_OPTION) $(COMM_OPTION)
### For compiling with acis 6.0 and higher, remove the -compat
### flag. (if you want compat add it after the SOLID_MODELER_DEFINES)
MACH_CXXFLAGS = $(ACIS_DEFINES) -DCUBIT -D__EXTERN_C__ -DSOLARIS \
                -DSUN -DSVR4 -DSYSV -DANSI -DXTFUNCPROTO +w \
		-DCUBIT_HOOPS_VERSION=$(CUBIT_HOOPS_VERSION) \
	        $(SOLID_MODELER_DEFINES) $(MEDIAL_DEFINE) \
                $(GRAPHICS_OPTION) $(COMM_OPTION)

# at least all -D and -I options from MACH_CXXFLAGS go here for makedepend
MACH_DEPEND_FLAGS = $(ACIS_DEFINES) -D__EXTERN_C__ -DSOLARIS \
                -DSUN -DSVR4 -DSYSV -DANSI -DXTFUNCPROTO +w \
		-DCUBIT_HOOPS_VERSION=$(CUBIT_HOOPS_VERSION) \
	        $(SOLID_MODELER_DEFINES) $(MEDIAL_DEFINE) \
                $(GRAPHICS_OPTION) $(COMM_OPTION)


# The -xpg is for profiling data, and if it isn't at the begining of this
# and you want to run gprof you'll need to add it to the begining of this 
# next line
MACH_LFLAGS = -R ${HOOPS_LIB_DIR}:/usr/openwin/lib

# Define the flags for building shared libraries
LD_SHARE = -G

# To compile with ACIS 6.0 and no compat, uncomment the following line and
# comment the one after.
#For ACIS 6.0 and later with no compat

#For a static link
# this next static link line doesn't work with cubit-beta and SUNWSpro.old
# static link causes problems accross operating systems
# dynamic appears better here.
#MACH_LIBS_LINK = -ldl -Bstatic -lnsl -lgen -lsocket\
#                 -lc -lCrun -lCstd -liostream -Bdynamic
#For a dynamic link
#MACH_LIBS_LINK =  -ldl -lnsl -Bstatic -lgen -Bdynamic -lsocket -lc \
#	-Bstatic -lCrun -lCstd -liostream  -Bdynamic
#Workshop 6 Update 1
#MACH_LIBS_LINK = -liostream
#Workshop 6 Update 2
MACH_LIBS_LINK = -liostream -lCstd

#For pre-ACIS 6.0 and compat
#MACH_LIBS_LINK = -KPIC -ldl -lnsl -lgen -lsocket\
#                 -lc -lCrun -lCstd -liostream
#-L$(SUN_WORKSHOP)/SC5.0/lib 
#MACH_LIBS_LINK = -ldl -lnsl -lgen -lsocket\
#                 -lc -lC
# MPI Communication stuff - just comment out if not making parallel version
#MPI_ARCH        = solaris
#MPI_COMM        = ch_p4
#MPI_INSTALL_DIR = /usr/netpub/ftp.mcs.anl.gov/pub/mpi/mpi-1.0.11-install
#MPI_LIB_PATH    = -L$(MPI_INSTALL_DIR)/lib/$(MPI_ARCH)/$(MPI_COMM)
#MPI_LIB_LIST    = -lmpi -lsocket -lnsl -lthread
#MPI_INCLUDE_DIR =  $(MPI_INSTALL_DIR)/include

# Some controls on echoing of compile commands; for verbose output, comment
# out the following two definitions
PREFIX = @
LINK_PREFIX = ${PREFIX}
ARCHIVE_PREFIX = ${PREFIX}
ECHO_COMMAND = @ echo "Compiling: $<"

## You can change this to filter the output of the compilation
#
#ACIS_GREP_FILTER = /usr/bin/egrep -v '(([.]hxx)|ATTRIB).*(::fixup_copy|::copy_common)'
#ACIS_FILTER = 2> compile_errors.txt ; \
#	${ACIS_GREP_FILTER} compile_errors.txt ; \
#	/bin/rm compile_errors.txt
ACIS_FILTER = 

# Some compiler definitions - check your environment to make sure the correct
# compilers are being referenced
#
# CC:
# for the SUN:
# (use the default, which is CC)
CC = cc
CXX = CC
#LINKER = purify ${CXX} -R$(ACIS_LIB_DIR)
#add compat to link line if ACIS < 6.0
LINKER = time ${CXX} -R$(ACIS_LIB_DIR)
ARCHIVER = ${CXX} -xar -o 

# MAKEDEENDEND:
# MAKEDEPEND:
# the following command uses the compiler to compute dependencies (uses option
# to exclude system includes)
#MAKEDEPEND = CC -xM1
# the following command uses our locally compiled makedepend to compute 
# dependencies (faster than compiler-driven make depend)
MAKEDEPEND = $(CUBITROOT)/bin/makedepend -f- -Y


CUBIT_HELP_DIR = $(SOURCE_DIR)
# you can change MAKE to make -j # for number of procs in your .cubit.myown file.
# MAKE = make -j 2
MAKE = make 
# Can use dmake on the Sun
#MAKE = dmake

#TEMPLATE_DIR = SunWS_cache

LEX = flex
YACC = bison
# YFLAGS is not needed with bison
#YFLAGS = -vdl

