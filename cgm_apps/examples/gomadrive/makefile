CUBIT_BASE_DIR = $(HOME)/cubit/clarolin/claro/cubit
CGM_BASE_DIR = $(CUBIT_BASE_DIR)
include $(CGM_BASE_DIR)/util/OtherFiles/make.client

CGM_GEOM_LIBS = $(GEOM_LIBS:%=$(CGM_BASE_DIR)/geom/%) \
                $(EXTRA_GEOM_LIBS:%=$(CGM_BASE_DIR)/geom/%)

CGM_LIBS_LINK = $(GEOM_LIBS_LINK:-L%=-L$(CGM_BASE_DIR)/%) \
                $(EXTRA_GEOM_LIBS_LINK:-Lgeom/%=-L$(CGM_BASE_DIR)/geom/%)

### main rule

drive: main.o gm_cgm_c_interface.o ${CGM_GEOM_LIBS}
	${LINKER} ${LFLAGS} ${MACH_LFLAGS} ${DEBUG_FLAG} main.o gm_cgm_c_interface.o ${CGM_LIBS_LINK} \
	-L${CUBIT_BASE_DIR}/geom/ACIS -lcubit_ACIS \
	-L${CUBIT_BASE_DIR}/geom/facet -lcubit_facet \
	-L${CUBIT_BASE_DIR}/geom/facetbool -lcubit_facetboolstub \
	-L${CUBIT_BASE_DIR}/geom -lcubit_geom \
	-L${CUBIT_BASE_DIR}/geom/Cholla -lcholla \
	-L${CUBIT_BASE_DIR}/util -lcubit_util \
	${CUBIT_SYS_LINK} -o drive 

$(CGM_GEOM_LIBS): 
	cd $(CGM_BASE_DIR); $(MAKE) $(@:$(CGM_BASE_DIR)/geom/%=%)

clean_all:
	find . -name "*.o" -exec rm -f \{\} \;
	rm -f drive

## replace some default suffix rules
.SUFFIXES : .y .o .l .cpp .cxx

.c.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CC} ${CFLAGS} ${MACH_CFLAGS} $(MOTIF_OPTIONS) \
        ${ALL_CUBIT_INCLUDES} -c $<

.cxx.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} -c $<

.cpp.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} -I. -DUSE_CGM -c $<

.C.o : 
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} $(MOTIF_OPTIONS) \
        ${CGM_INCLUDES} -I. -DUSE_CGM -c $<

.l.c:
	${ECHO_COMMAND}
	${LEX} -t ${L_FLAGS} $< > $@

.y.o:
	${ECHO_COMMAND}
	${YACC} ${YFLAGS} $<
	${ECHO_COMMAND}
	${PREFIX} ${CXX} ${CXXFLAGS} ${MACH_CXXFLAGS} ${ALL_CUBIT_INCLUDES} -c y.tab.c
	mv y.tab.o $@



PREFIX =
