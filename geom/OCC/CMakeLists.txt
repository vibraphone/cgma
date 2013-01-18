
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

ADD_DEFINITIONS(
  -DHAVE_OCC_IGES
  -DHAVE_OCC_STEP
  )

if(WIN32)
  ADD_DEFINITIONS(
  -DWNT
  -DHAVE_OCC_DEF
  -DOCC_INC_FLAG
  -DWNT
  )
endif()

# build cubit occ library from sources
ADD_LIBRARY(cubit_occ MODULE ${OCC_SRCS})
TARGET_LINK_LIBRARIES(cubit_occ cubit_geom cubit_util TKShHealing TKBinL TKBin TKLCAF TKCAF TKCDF TKPShape TKFillet TKOffset TKMesh TKBool TKBO TKFeat TKG2d TKG3d TKXSBase TKTopAlgo TKPrim TKGeomAlgo TKBRep TKGeomBase TKernel TKMath TKAdvTools TKIGES TKSTEP)


INSTALL(TARGETS cubit_occ ${CUBIT_GEOM_EXPORT_GROUP}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
        )
