PROJECT(cubit_occ)

if(CMAKE_SYSTEM_NAME MATCHES Windows)
  if(CMAKE_SIZEOF_VOID_P MATCHES 8)
    set(OCE_MODULE_HINTS "${CUBITROOT}/OCE/OCE-0.10-64/cmake")
  else()
    set(OCE_MODULE_HINTS "${CUBITROOT}/OCE/OCE-0.10-32/cmake")
  endif()
endif()

if(CMAKE_SYSTEM_NAME MATCHES Linux)
  if(CMAKE_SIZEOF_VOID_P MATCHES 8)
    set(OCE_MODULE_HINTS "${CUBITROOT}/OCE/OCE-0.10-64/lib/oce-0.10")
  else()
    set(OCE_MODULE_HINTS "${CUBITROOT}/OCE/OCE-0.10-32/lib/oce-0.10")
  endif()
endif()

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(OCE_MODULE_HINTS "${CUBITROOT}/OCE/OCE-0.10/OCE.framework/Versions/0.10/Resources")
endif()

find_package(OCE CONFIG HINTS ${OCE_MODULE_HINTS})

include_directories(${OCE_INCLUDE_DIRS})
set(OCE_INCLUDE_DIRS ${OCE_INCLUDE_DIRS} PARENT_SCOPE)

SET(OCC_SRCS
  OCCAttribSet.cpp
  OCCAttribSet.hpp
  OCCBody.cpp
  OCCBody.hpp
  OCCCoEdge.cpp
  OCCCoEdge.hpp
  OCCCoFace.cpp
  OCCCoFace.hpp
  OCCCurve.cpp
  OCCCurve.hpp
  OCCDrawTool.cpp
  OCCDrawTool.hpp
  OCCGeometryCreator.hpp
  OCCHistory.cpp
  OCCHistory.hpp
  OCCLoop.cpp
  OCCLoop.hpp
  OCCLump.cpp
  OCCLump.hpp
  OCCModifyEngine.cpp
  OCCModifyEngine.hpp
  OCCPoint.cpp
  OCCPoint.hpp
  OCCQueryEngine.cpp
  OCCQueryEngine.hpp
  OCCShapeAttributeSet.cpp
  OCCShapeAttributeSet.hpp
  OCCShell.cpp
  OCCShell.hpp
  OCCSurface.cpp
  OCCSurface.hpp
)

set(OCC_HEADERS
  OCCAttribSet.hpp
  OCCBody.hpp
  OCCCoEdge.hpp
  OCCCoFace.hpp
  OCCCurve.hpp
  OCCDrawTool.hpp
  OCCGeometryCreator.hpp
  OCCLoop.hpp
  OCCLump.hpp
  OCCModifyEngine.hpp
  OCCPoint.hpp
  OCCQueryEngine.hpp
  OCCShapeAttributeSet.hpp
  OCCShell.hpp
  OCCSurface.hpp
)

# FIXME: Detect presence of OCC libraries for IGES and STEP
set(CGMA_OCC_DEFINES "-DHAVE_OCC -DHAVE_OCC_IGES -DHAVE_OCC_STEP -DHAVE_OCC_STL" PARENT_SCOPE)
add_definitions(${CGMA_OCC_DEFINES})

if(WIN32)
  ADD_DEFINITIONS(
  -DWNT
  -DHAVE_OCC_DEF
  -DOCC_INC_FLAG
  -DWNT
  )
endif()

include(${OCE_DIR}/OCEConfig.cmake)
include(${OCE_DIR}/OCE-libraries.cmake)
# build cubit occ library from sources
set(CGMA_OCC_DEPLIBS
  #TKShHealing TKBinL TKBin TKLCAF TKCAF TKCDF
  #TKPShape TKFillet TKOffset TKMesh TKBool TKBO TKFeat
  #TKG2d TKG3d TKXSBase TKTopAlgo TKPrim TKGeomAlgo TKBRep TKGeomBase
  #TKernel TKMath TKAdvTools TKIGES TKSTEP TKSTL

  #TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKMesh TKPrim TKBO TKBool TKFeat TKFillet
  #TKShHealing TKOffset TKHLR TKCAF TKCDF TKLCAF TKBinL TKXSBase TKIGES TKXSBase TKSTEPBase TKSTEPAttr TKSTEP209
  #TKSTEP TKSTL

  #${OCE_LIBRARIES}

  PARENT_SCOPE
)
add_library(cubit_occ SHARED ${OCC_SRCS})
target_link_libraries(cubit_occ
  cubit_geom
  cubit_util
  ${OCE_LIBRARIES}
)
set(CGMA_LIBS "${CGMA_LIBS}" cubit_occ)
set(CGMA_LIBS "${CGMA_LIBS}" PARENT_SCOPE)

install(
  TARGETS cubit_occ ${CUBIT_GEOM_EXPORT_GROUP}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
)
install(
  FILES ${OCC_HEADERS}
  DESTINATION include
  COMPONENT Development
)
