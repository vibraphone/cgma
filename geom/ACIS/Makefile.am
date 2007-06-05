# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFS = $(CANT_USE_STD) $(CANT_USE_STD_IO) $(TEMPLATE_DEFS_INCLUDED) \
       $(LITTLE_ENDIAN) $(ACIS_STEP_TRANSLATOR) $(ACIS_IGES_TRANSLATOR) \
       -DCUBIT_ACIS_VERSION=$(ACIS_VERSION) -DACIS_VERSION=$(ACIS_VERSION) \
       -D$(ACIS_PLATFORM) -DACIS_LOCAL_OPS
LIBS := -L$(ACIS_LIB_DIR) $(ACIS_LIBS) $(LIBS)
INCLUDES = -I$(top_srcdir)/util \
           -I$(top_srcdir)/geom \
           -I$(ACIS_DIR)/include

# The name of the library to build
lib_LTLIBRARIES = 
if BUILD_CGM
  lib_LTLIBRARIES += libcubit_ACIS.la
endif

# The directory where headers will be installed
libcubit_ACIS_la_includedir = $(includedir)

# The non-template sources
libcubit_ACIS_la_SOURCES = \
    AcisBridge.cpp \
    AcisEdgeTool.cpp \
    AcisFacetManager.cpp \
    AcisFeatureEngine.cpp \
    AcisHealerTool.cpp \
    AcisModifyEngine.cpp \
    AcisQueryEngine.cpp \
    AcisSurfaceTool.cpp \
    AcisToolUtil.cpp \
    AcisTweakTool.cpp \
    AcisDrawTool.cpp \
    attrib_cubit_owner.cpp \
    attrib_gtc.cpp \
    attrib_gtc_name.cpp \
    attrib_snl.cpp \
    attrib_snl_simple.cpp \
    BodyACIS.cpp \
    CoEdgeACIS.cpp \
    CurveACIS.cpp \
    LoopACIS.cpp \
    LumpACIS.cpp \
    PointACIS.cpp \
    ShellACIS.cpp \
    SurfaceACIS.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_ACIS_la_include_HEADERS = \
    AcisBridge.hpp \
    AcisDrawTool.hpp \
    AcisEdgeTool.hpp \
    AcisFacetManager.hpp \
    AcisFeatureEngine.hpp \
    AcisHealerTool.hpp \
    AcisModifyEngine.hpp \
    AcisQueryEngine.hpp \
    AcisSurfaceTool.hpp \
    AcisToolUtil.hpp \
    AcisTweakTool.hpp \
    AcisTypes.h \
    BodyACIS.hpp \
    CoEdgeACIS.hpp \
    CurveACIS.hpp \
    LoopACIS.hpp \
    LumpACIS.hpp \
    PointACIS.hpp \
    ShellACIS.hpp \
    SurfaceACIS.hpp \
    attrib_cubit_owner.hpp \
    attrib_gtc.hpp \
    attrib_gtc_name.hpp \
    attrib_snl.hpp \
    attrib_snl_simple.hpp \
    decl_none.h
  
