CGM_DIR = ../../..
include ${CGM_DIR}/lib/cgm.make

all: webcut mergechk imprint

webcut: webcut.o  ${CGM_LIBS_FILES}
	$(CGM_CXX) -o webcut $(LDFLAGS) webcut.o $(CGM_LIBS_LINK) 

mergechk:  mergechk.o  ${CGM_LIBS_FILES}
	$(CGM_CXX) -o mergechk $(LDFLAGS) mergechk.o $(CGM_LIBS_LINK)

imprint:  imprint.o ${CGM_LIBS_FILES}
	$(CGM_CXX) -o  imprint $(LDFLAGS) imprint.o $(CGM_LIBS_LINK)

.cpp.o:
	$(CGM_CXX) ${CGM_CXXFLAGS} ${CGM_INCLUDES} -c -o $@ $<

