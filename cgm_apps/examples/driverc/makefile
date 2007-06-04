CUBIT_BASE_DIR = /export/home/tjtautg/cubit/work
CGM_BASE_DIR = /export/home/tjtautg/cubit/cgm
include ${CUBIT_BASE_DIR}/util/OtherFiles/make.client

<<<<<<< makefile
mergechkc: main.o mergechkc.o $(CUBIT_LIBS) ${CCAPI_LIBS}
	${LINKER} ${DEBUG_FLAG} ${LFLAGS} main.o mergechkc.o ${CCAPI_LIBS_LINK} ${GEOM_LIBS_LINK} \
                   ${ACIS_LINK} ${CUBIT_SYS_LINK} -o mergechkc

=======
mergechkc: mergechkc.o $(CUBIT_LIBS) ${CCAPI_LIBS}
	${LINKER} ${DEBUG_FLAG} ${LFLAGS} mergechkc.o ${CCAPI_LIBS_LINK} ${GEOM_LIBS_LINK} \
                   ${ACIS_LINK} ${CUBIT_SYS_LINK} -o mergechkc

>>>>>>> 1.2
driverc: main.o driverc.o $(CUBIT_LIBS) ${CCAPI_LIBS}
	${LINKER} ${DEBUG_FLAG} ${LFLAGS} main.o driverc.o ${CCAPI_LIBS_LINK} ${GEOM_LIBS_LINK} \
                   ${ACIS_LINK} ${CUBIT_SYS_LINK} -o driverc

libs: 
	cd $(CUBIT_BASE_DIR); make ${GEOM_LIBS}
	cd $(CUBIT_BASE_DIR); make ${CCAPI_LIBS}

apilibs: ${GEOM_LIBS}
	$(MAKE) libs
	cd $(CUBIT_BASE_DIR)/util/ccapi; make
	cd $(CUBIT_BASE_DIR)/geom/ccapi; make

depend : driverc.c mergechkc.c
	${MAKEDEPEND} -f make.dependencies \
           ${MACH_CXXFLAGS} ${MACH_CFLAGS}\
           ${CCAPI_INCLUDES} ${ALL_CUBIT_INCLUDES} driverc.c mergechkc.c ${MAKE_DEPEND_INCLUDE} 

## replace some default suffix rules
.SUFFIXES : .y .o .l .cpp .cxx .c 

.c.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CC} ${DEBUG_FLAG} ${CFLAGS} ${MACH_CFLAGS} $(MOTIF_OPTIONS) \
        ${CCAPI_INCLUDES} ${ALL_CUBIT_INCLUDES} -c $<

.cxx.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${DEBUG_FLAG} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${ALL_CUBIT_INCLUDES} -c $<

.cpp.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${DEBUG_FLAG} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${ALL_CUBIT_INCLUDES} -c $<

.l.c:
	${ECHO_COMMAND}
	${LEX} -t ${L_FLAGS} $< > $@

.y.o:
	${ECHO_COMMAND}
	${YACC} ${YFLAGS} $<
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} ${ALL_CUBIT_INCLUDES} -c y.tab.c
	mv y.tab.o $@


