PROJECT(cubit_occ)

SET(OCC_SRCS
    OCCShapeAttributeSet.cpp
    OCCAttribSet.cpp
    OCCBody.cpp
    OCCCoEdge.cpp
    OCCCoFace.cpp
    OCCCurve.cpp
    OCCLoop.cpp
    OCCLump.cpp
    OCCModifyEngine.cpp
    OCCPoint.cpp
    OCCQueryEngine.cpp
    OCCShell.cpp
    OCCSurface.cpp
    OCCDrawTool.cpp
   )

SET(OCC_HDRS
    OCCShapeAttributeSet.hpp
    OCCAttribSet.hpp
    OCCBody.hpp
    OCCCoEdge.hpp
    OCCCoFace.hpp
    OCCCurve.hpp
    OCCGeometryCreator.hpp
    OCCLoop.hpp
    OCCLump.hpp
    OCCModifyEngine.hpp
    OCCPoint.hpp
    OCCQueryEngine.hpp
    OCCShell.hpp
    OCCSurface.hpp
    OCCDrawTool.hpp
    )

SET(OCC_LIBS
     TKMath
     TKXSBase
     TKernel
     TKShHealing
     TKBRep
     TKG3d
     TKBO
     TKFeat
     TKFillet
     TKG2d
     TKOffset
     TKBool
     TKGeomAlgo
     TKGeomBase
     TKPrim
     TKTopAlgo
     TKHLR
     TKIGES
     TKMesh
     TKPrim
     TKSTEP209
     TKSTEPAttr
     TKSTEPBase
     TKSTEP
     TKSTL
     TKTopAlgo
     TKXSBase
     TKLCAF
     TKBinL
    )

SET(OCC_DIR "" CACHE PATH "Location where OpenCASCADE in installed")
FIND_LIBRARY( OCC_LIBRARY TKernel ${OCC_DIR}/lib )
FIND_PATH( OCC_INCLUDES Standard_Version.hxx ${OCC_DIR}/inc )

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif()
ADD_LIBRARY( cubit_OCC ${OCC_SRCS} ${OCC_HDRS} )
INCLUDE_DIRECTORIES( ${OCC_INCLUDES} )
TARGET_LINK_LIBRARIES( cubit_OCC ${OCC_LIBRARY} ${OCC_LIBS} )
SET_SOURCE_FILES_PROPERTIES( ${OCC_SRCS} COMPILE_FLAGS "-DHAVE_IOSTREAM -DHAVE_IOMANIP -DHAVE_FSTREAM -DHAVE_LIMITS_H" )
