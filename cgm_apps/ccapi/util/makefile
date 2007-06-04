CUBIT_BASE_DIR = /export/home/tjtautg/cubit/work
CGM_BASE_DIR = /export/home/tjtautg/cubit/cgm
include ${CUBIT_BASE_DIR}/util/OtherFiles/make.client

SHELL = /bin/sh

UTIL_CCAPI_SRCS = \
	CubitApp_ccapi.cpp \
	CubitObserver_ccapi.cpp \
	CubitUtil_ccapi.cpp

ALL_HPPS = ${UTIL_CCAPI_SRCS:.cpp=.h}

### blank rule to force rebuilds of the main cubit target
default: libcubit_util_ccapi.a 

force_rebuild :

libcubit_util_ccapi.a : ${UTIL_CCAPI_SRCS:.cpp=.o} makefile
	${ARCHIVER} libcubit_util_ccapi.a ${UTIL_CCAPI_SRCS:.cpp=.o}

### utility targets

clean_all :
	@- rm -f *.o *~ #* 
	@- rm -rf ${TEMPLATE_DIR}
	@- rm -f ${SOURCE_DIR}/*~ 
	@- rm -f ${SOURCE_DIR}/#*

clean :
	@- rm -f  *~ #* 
	@- rm -f ${SOURCE_DIR}/*~ 
	@- rm -f ${SOURCE_DIR}/#*

## top target to make dependencies
depend : ${UTIL_CCAPI_SRCS}
	${MAKEDEPEND} -f make.dependencies \
           ${MACH_CXXFLAGS} ${MACH_CFLAGS}\
           ${CGM_INCLUDES} ${CCAPI_INCLUDES} ${CUBIT_SYS_INCLUDES} ${UTIL_CCAPI_SRCS} ${MAKE_DEPEND_INCLUDE} 

## replace some default suffix rules
.SUFFIXES : .y .o .l .cpp .cxx

.c.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CC} ${DEBUG_FLAG} ${CFLAGS} ${MACH_CFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} ${CCAPI_INCLUDES} ${CUBIT_SYS_INCLUDES} -c $<

.cxx.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${DEBUG_FLAG} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} ${CCAPI_INCLUDES} ${CUBIT_SYS_INCLUDES} -c $<

.cpp.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${DEBUG_FLAG} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} ${CCAPI_INCLUDES} ${CUBIT_SYS_INCLUDES} -c $<

.l.c:
	${ECHO_COMMAND}
	${LEX} -t ${L_FLAGS} $< > $@

.y.o:
	${ECHO_COMMAND}
	${YACC} ${YFLAGS} $<
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} ${ALL_INCLUDES} -c y.tab.c
	mv y.tab.o $@

include make.dependencies

# DO NOT DELETE THIS LINE -- make depend depends on it.


