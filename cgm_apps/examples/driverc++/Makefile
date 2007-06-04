CGM_DIR = ../../..
include ${CGM_DIR}/config.make

all: mergechk

mergechk: mergechk.o ${CGM_LIBS_FILES}
	$(CXX) $(LDFLAGS) mergechk.o $(CGM_LIBS_LINK)

.cpp.o:
	$(CXX) ${CGM_CXXFLAGS} ${CGM_INCLUDES} -c -o $@ $<

