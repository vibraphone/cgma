
PROJECT(cubit_facetbool)

set(FACETBOOL_SRCS
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

set(FACETBOOL_HEADERS
  FBDefines.hpp
  FBStructs.hpp
)
foreach(var ${FACETBOOL_SRCS})
  string(REGEX REPLACE ".cpp" ".hpp" header ${var})
  set(FACETBOOL_HEADERS ${FACETBOOL_HEADERS} ${header})
endforeach(var ${FACETBOOL_SRCS})

add_library(cubit_facetboolstub STATIC ${FACETBOOL_SRCS})
set(CGMA_LIBS "${CGMA_LIBS}" cubit_facetboolstub PARENT_SCOPE)

install(
  TARGETS cubit_facetboolstub
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)
install(
  FILES ${FACETBOOL_HEADERS}
  DESTINATION include
  COMPONENT Development
)
