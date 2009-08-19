#
# CGM_BASE_DIR: base directory below which are geom/ and util/, the cgm source
#
CGM_BASE_DIR = ${HOME}/lib/CGM
#
# CUBIT_INSTALL_DIR: the installation directory of your cubit executable, below
# which should be bin/ containing libcubiti19.so (in bin/) and other 
# shared libs (in libs/).
#
CUBIT_INSTALL_DIR = ${HOME}/cubit-10.2-amd64

LFLAGS = -Wl,--unresolved-symbols=ignore-all -Wl,--allow-shlib-undefined -Wl,-rpath,${CUBIT_INSTALL_DIR}/bin -Wl,-rpath,${HOME}/lib/OpenCASCADE6.3.0/lib

include ${CGM_BASE_DIR}/lib/cgm.make

cgm_cubit_acis_app: main_cubit_acis.o Makefile snippet.o
	${CXX} ${CGM_INCLUDES} -o cgm_cubit_acis ${LFLAGS} main_cubit_acis.o snippet.o ${CGM_LIBS_LINK}

main_cubit_acis.o : Makefile NoAcisQueryEngine.hpp NoAcisModifyEngine.hpp

mcnp_converter.o : mcnp_converter.cpp Makefile NoAcisQueryEngine.hpp NoAcisModifyEngine.hpp

main_acis.o: main_acis.cpp
	$(CXX) ${CGM_CXXFLAGS} -g -I. ${CGM_INCLUDES} -o $@ -c $<

test:
	./cgm_cubit_acis ${HOME}/source/CGM/itaps/testgeom.sat

.cpp.o:
	$(CXX) ${CGM_CXXFLAGS} -g -I. ${CGM_INCLUDES} -o $@ -c $<

