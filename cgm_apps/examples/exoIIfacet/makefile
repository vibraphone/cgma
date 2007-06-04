CUBIT_BASE_DIR = $(HOME)/cubit/cgm
CGM_BASE_DIR = $(CUBIT_BASE_DIR)
include $(CUBIT_BASE_DIR)/util/OtherFiles/make.client

CGM_GEOM_LIBS = $(GEOM_LIBS:%=$(CGM_BASE_DIR)/geom/%) \
                $(EXTRA_GEOM_LIBS:%=$(CGM_BASE_DIR)/geom/%) 

CGM_LIBS_LINK = $(GEOM_LIBS_LINK:-L%=-L$(CGM_BASE_DIR)/%) \
                $(EXTRA_GEOM_LIBS_LINK:-Lgeom/%=-L$(CGM_BASE_DIR)/geom/%) 

### main rule

exo: ExoLoader.o ElemConnectivity.o $(CGM_GEOM_LIBS)
	${LINKER} ${LFLAGS} ${MACH_LFLAGS} ${DEBUG_FLAG} \
	ExoLoader.o ElemConnectivity.o \
	${CGM_LIBS_LINK} \
	-L${CUBIT_BASE_DIR}/geom/facetbool -lcubit_facetboolstub \
	${CUBIT_SYS_LINK} \
	-L${EXODUS_LIB_DIR} -lexoIIv2c \
	-L$(NETCDF_LIB_DIR) -lnetcdf -o exo 

$(CGM_GEOM_LIBS): 
	cd $(CGM_BASE_DIR); $(MAKE) $(@:$(CGM_BASE_DIR)/geom/%=%)

clean_all:
	find . -name "*.o" -exec rm -f \{\} \;
	rm -f exo

## replace some default suffix rules
.SUFFIXES : .y .o .l .cpp .cxx

.c.o : 
	${CC} ${CFLAGS} ${MACH_CFLAGS} $(MOTIF_OPTIONS) \
        ${ALL_CUBIT_INCLUDES} -c $<

.cxx.o : 
	${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} ${EXODUS_INCLUDE} ${NETCDF_INCLUDE} -c $<

.cpp.o : 
	${CXX} -g ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} ${EXODUS_INCLUDE} ${NETCDF_INCLUDE} -c $<

.l.c:
	${ECHO_COMMAND}
	${LEX} -t ${L_FLAGS} $< > $@

.y.o:
	${ECHO_COMMAND}
	${YACC} ${YFLAGS} $<
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} ${ALL_CUBIT_INCLUDES} -c y.tab.c
	mv y.tab.o $@



