
PROJECT(cubit_acis)

SET(ACIS_SRCS
    AcisBridge.cpp
    AcisEdgeTool.cpp
    AcisFacetManager.cpp
    AcisFeatureEngine.cpp
    AcisHealerTool.cpp
    AcisModifyEngine.cpp
    AcisQueryEngine.cpp
    AcisSurfaceTool.cpp
    AcisToolUtil.cpp
    AcisTweakTool.cpp
    AcisDrawTool.cpp
    attrib_cubit_owner.cpp
    attrib_snl.cpp
    attrib_snl_simple.cpp
    BodyACIS.cpp
    CoEdgeACIS.cpp
    CurveACIS.cpp
    LoopACIS.cpp
    LumpACIS.cpp
    PointACIS.cpp
    ShellACIS.cpp
    SurfaceACIS.cpp
    )

# for convenience, include this in the project
IF(WIN32)
  FOREACH(var ${ACIS_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(ACIS_HEADERS ${ACIS_HEADERS} ${header})
  ENDFOREACH(var ${ACIS_SRCS})
  IF(NOT ${CMAKE_GENERATOR} MATCHES "NMake")
    SET(EXTRA_ACIS_SRCS UseACIS.cmake)
  ENDIF(NOT ${CMAKE_GENERATOR} MATCHES "NMake")
ENDIF(WIN32)


INCLUDE(${cubit_acis_SOURCE_DIR}/UseACIS.cmake)

# add acis includes to current project
INCLUDE_DIRECTORIES(${ACIS_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/gtcAttrib/incl)

# build cubit acis library from sources
ADD_LIBRARY(cubit_ACIS ${ACIS_SRCS} ${ACIS_HEADERS} ${EXTRA_ACIS_SRCS})

SUBDIRS(gtcAttrib)

IF(UNIX)
ADD_CUSTOM_TARGET(switch_acis ${CMAKE_MAKE_PROGRAM} clean)
ENDIF(UNIX)

# install ACIS libs 
SET(CMAKE_INSTALL_3RD_PARTY_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "Install prefix for 3rd party libraries")
SET(CMAKE_INSTALL_BINARY_DIR "bin" CACHE PATH "Install directory for binaries")
IF(CMAKE_INSTALL_3RD_PARTY_PREFIX)
  SET(ACIS_INSTALL_PREFIX ${CMAKE_INSTALL_3RD_PARTY_PREFIX}/${CMAKE_INSTALL_BINARY_DIR})
  IF(WIN32)
    STRING(REGEX REPLACE "/lib/" "/bin/" ACIS_SHARED_LIBRARY_DIR ${ACIS_LIBRARY_DIR})
    STRING(REGEX REPLACE "/lib/" "/bin/" ACIS_SHARED_LIBRARY_DIR_DEBUG ${ACIS_LIBRARY_DIR_DEBUG})
    SET(ACIS_INSTALL_LIBS ${ACIS_RELEASE_LIBS})
    SET(ACIS_INSTALL_LIBS_DEBUG ${ACIS_DEBUG_LIBS})
  ELSE(WIN32)
    SET(ACIS_SHARED_LIBRARY_DIR ${ACIS_LIBRARY_DIR})
    SET(ACIS_INSTALL_LIBS ${ACIS_LIBS})
  ENDIF(WIN32)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/InstallACIS.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/InstallACIS.cmake @ONLY)
  SET_TARGET_PROPERTIES(cubit_ACIS PROPERTIES POST_INSTALL_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/InstallACIS.cmake)
ENDIF(CMAKE_INSTALL_3RD_PARTY_PREFIX)

