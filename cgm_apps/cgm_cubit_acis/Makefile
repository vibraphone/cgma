#
# CGM_BASE_DIR: base directory below which are geom/ and util/, the cgm source
#
CGM_BASE_DIR = ${HOME}/cubit/CGM-10.1
#
# CUBIT_INSTALL_DIR: the installation directory of your cubit executable, below
# which should be bin/ and libs/ containing libcubiti19.so (in bin/) and the acis 
# shared libs (in libs/).
#
CUBIT_INSTALL_DIR = /cubit

LFLAGS = -Wl,--unresolved-symbols=ignore-all -Wl,--allow-shlib-undefined -Wl,-rpath,${CUBIT_INSTALL_DIR}/bin -Wl,-rpath,${CUBIT_INSTALL_DIR}/libs

cgm_cubit_acis_app: main_cubit_acis.o Makefile snippet.o
	${CXX} -o cgm_cubit_acis ${LFLAGS} main_cubit_acis.o snippet.o -L${CUBIT_INSTALL_DIR}/bin -L${CUBIT_INSTALL_DIR}/libs \
	-lcubiti19 -lacisstep -lxstep -lacisiges -lxiges -lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXBase -lxacis2k -lxcore2k -lSpaAVis -lSpaAWarp -lSpaASurf -lSpaALops -lSpaABlend -lSpaACIS -lSpaBase -lcubiti19 -lacisstep -lxstep -lacisiges -lxiges -lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXBase -lxacis2k -lxcore2k -lSpaAVis -lSpaAWarp -lSpaASurf -lSpaALops -lSpaABlend -lSpaACIS -lSpaBase

include ${CGM_BASE_DIR}/include/config.make

main_cubit_acis.o : Makefile NoAcisQueryEngine.hpp NoAcisModifyEngine.hpp

mcnp_converter.o : mcnp_converter.cpp Makefile NoAcisQueryEngine.hpp NoAcisModifyEngine.hpp

main_acis.o: main_acis.cpp
	$(CXX) ${CGM_CXXFLAGS} -g -I. ${CGM_INCLUDES} -o $@ -c $<

test:
	./cgm_cubit_acis geom.sat

.cpp.o:
	$(CXX) ${CGM_CXXFLAGS} -g -I. -I${CGM_BASE_DIR}/include -o $@ -c $<

