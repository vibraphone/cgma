
PROJECT(cubit_facetbool)

SET(FACETBOOL_SRCS
    FBClassify.cpp
    FBDataUtil.cpp
    FBImprint.cpp
    FBIntersect.cpp
    FBPolyhedron.cpp
    FBRetriangulate.cpp
    FBTiler.cpp
    IntegerHash.cpp
    KdTree.cpp
    )
    
IF(WIN32)
  FOREACH(var ${FACETBOOL_SRCS})
    STRING(REGEX REPLACE ".cpp" ".hpp" header ${var})
    SET(FACETBOOL_HEADERS ${FACETBOOL_HEADERS} ${header})
  ENDFOREACH(var ${FACETBOOL_SRCS})
ENDIF(WIN32)

ADD_LIBRARY(cubit_facetboolstub ${FACETBOOL_SRCS} ${FACETBOOL_HEADERS})
