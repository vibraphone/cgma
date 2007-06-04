CGM_DIR = ../../..
include $(CGM_DIR)/config.make

### main rule

facets: facets.o ${CGM_LIBS_FILES}
	${CXX} ${LDFLAGS} facets.o ${CGM_LIBS_LINK} -o facets

check: facets

check-TESTS: check
	./facets 90

clean_all:
	find . -name "*.o" -exec rm -f \{\} \;
	rm -f facets

## replace some default suffix rules
.SUFFIXES : .o .cpp

.cpp.o : 
	${CXX} ${CXXFLAGS} ${CGM_INCLUDES} -c $<
