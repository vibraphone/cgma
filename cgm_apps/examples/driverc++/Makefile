CGM_DIR = ../../..
include ${CGM_DIR}/lib/cgm.make

all: mergechk 

mergechk:  mergechk.o  ${CGM_LIBS_FILES}
	$(CGM_CXX) -o mergechk $(LDFLAGS) mergechk.o $(CGM_LIBS_LINK)

.cpp.o:
	$(CGM_CXX) ${CGM_CXXFLAGS} ${CGM_INCLUDES} -c -o $@ $<

