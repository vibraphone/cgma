option(ITAPS_USE_BABEL "Should ITAPS use Babel to generate C and FORTRAN clients?" "OFF")

if (ITAPS_USE_BABEL)
  add_subdirectory(SIDL)
endif()

set(CGMA_IGEOM_SRCS
	CATag.cpp
	iGeom_CGMA.cc
	iGeomError.cc
)

set(CGMA_IGEOM_HDRS
  iBase.h
  iBase_f.h
  iGeom.h
  iGeom_FCDefs.h
  iGeom_f.h
  #iGeom_protos.h
)

include_directories(
  ${cubit_geom_SOURCE_DIR}
  ${cubit_util_SOURCE_DIR}
  ${cubit_geom_SOURCE_DIR}/ACIS
  ${cubit_geom_SOURCE_DIR}/virtual
  ${cubit_geom_SOURCE_DIR}/facet
  ${cubit_geom_SOURCE_DIR}/cholla
  ${cgma_init_SOURCE_DIR}
  ${cubit_geom_BINARY_DIR}
)

add_library(iGeom ${CGMA_IGEOM_SRCS})
target_link_libraries(iGeom 
  LINK_PUBLIC
    cgm
)

install(
  TARGETS iGeom
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)
install(
  FILES ${CGMA_IGEOM_HDRS}
  DESTINATION include
  COMPONENT Development
)
