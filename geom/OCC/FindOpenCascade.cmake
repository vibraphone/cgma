#
# ======================================================================================
# OpenCascade related stuff
# sets OPENCASCADE_INCLUDE_DIRS and OPENCASCADE_LIBRARIES
# ======================================================================================


find_path(OPENCASCADE_INC_DIR NAMES PGeom_Circle.hxx PATH_SUFFIXES opencascade)
#mark_as_advanced(OPENCASCADE_INC_DIR)
set(OPENCASCADE_INCLUDE_DIRS ${OPENCASCADE_INC_DIR})
set(OPENCASCADE_LIBRARIES)

set(OPENCASCADE_FOUNDATION_LIBS TKernel TKMath TKAdvTools TKjcas)
set(OPENCASCADE_MODELING_LIBS  TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo
                              TKTopAlgo TKPrim TKBO TKHLR TKMesh TKShHealing
                              TKBool TKXMesh TKFillet TKFeat TKOffset
                              TKSTL TKXSBase TKSTEPBase TKIGES TKSTEPAttr
                              TKSTEP209 TKSTEP TKCDF TKPShape TKLCAF TKCAF TKBinL TKBin)
set(OPENCASCADE_OCAF_LIBS  TKLCAF TKBin TKXml TKBinTObj TKXmlTObj TKPCAF
                          TKStdSchema StdPlugin XmlPlugin BinPlugin BinTObjPlugin
                          XmlTObjPlugin TKXCAF TKXCAFSchema TKXmlXCAF TKBinXCAF
                          TKXDEIGES TKXDESTEP XCAFPlugin XmlXCAFPlugin
                          BinXCAFPlugin TKXmlL)
set(OPENCASCADE_VIS_LIBS TKService TKV2d TKV3d TKOpenGl TKMeshVS TKNIS TKVRML)


set(OPENCASCADE_LIB_DIR "" CACHE PATH "")
foreach(lib ${OPENCASCADE_MODELING_LIBS} ${OPENCASCADE_FOUNDATION_LIBS}) # ${OPENCASCADE_OCAF_LIBS})
  find_library(OPENCASCADE_${lib}_LIBRARY NAMES ${lib} PATHS ${OPENCASCADE_LIB_DIR})
  mark_as_advanced(OPENCASCADE_${lib}_LIBRARY)
  add_library(${lib} UNKNOWN IMPORTED)
  set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION "${OPENCASCADE_${lib}_LIBRARY}")
  set(OPENCASCADE_LIBRARIES ${OPENCASCADE_LIBRARIES} ${lib})
endforeach(lib)

   #message(STATUS "OpenCascade lib dir = ${OPENCASCADE_LIB_DIR}")
   #message(STATUS "OpenCascade inc dir = ${OPENCASCADE_INC_DIR}")
   #message(STATUS "OpenCascade libs    = ${OPENCASCADE_LIBRARIES}")


