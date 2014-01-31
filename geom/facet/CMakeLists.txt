project(cubit_facet)

set(FACET_SRCS
  FacetAttrib.cpp
  FacetAttribSet.cpp
  FacetBody.cpp
  FacetboolInterface.cpp
  FacetCoEdge.cpp
  FacetCurve.cpp
  #FacetGeometryCreator.cpp
  FacetLoop.cpp
  FacetLump.cpp
  FacetModifyEngine.cpp
  FacetParamTool.cpp
  FacetPoint.cpp
  FacetQueryEngine.cpp
  FacetShell.cpp
  FacetSurface.cpp
  GridSearchTree.cpp
)

set(FACET_HEADERS
  FacetAttrib.hpp
  FacetAttribSet.hpp
  FacetBody.hpp
  FacetCoEdge.hpp
  FacetCurve.hpp
  FacetGeometryCreator.hpp
  FacetLoop.hpp
  FacetLump.hpp
  FacetModifyEngine.hpp
  FacetParamTool.hpp
  FacetPoint.hpp
  FacetQueryEngine.hpp
  FacetShell.hpp
  FacetSurface.hpp
  FacetboolInterface.hpp
  GridSearchTree.hpp
  GridSearchTreeNode.hpp
)

include_directories(
  ${cubit_geom_SOURCE_DIR}/facetbool
  ${cubit_geom_SOURCE_DIR}/Cholla
)

add_library(cubit_facet STATIC ${FACET_SRCS} ${FACET_HEADERS})
set(CGMA_LIBS "${CGMA_LIBS}" cubit_facet PARENT_SCOPE)

install(
  TARGETS cubit_facet
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)
install(
  FILES ${FACET_HEADERS}
  DESTINATION include
  COMPONENT Development
)
