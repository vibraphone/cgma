project( cgma_init )

include_directories(
  ${cubit_geom_SOURCE_DIR}
  ${cubit_util_SOURCE_DIR}
  ${cubit_geom_SOURCE_DIR}/ACIS
  ${cubit_geom_SOURCE_DIR}/virtual
  ${cubit_geom_SOURCE_DIR}/facet
  ${cubit_geom_SOURCE_DIR}/cholla
  ${cubit_geom_SOURCE_DIR}/OCC
  ${cgma_init_SOURCE_DIR}
  ${cubit_geom_BINARY_DIR}
  ${OCE_INCLUDE_DIRS}
)

set(CGMA_INIT_SRCS
  InitCGMA.cpp
)
set(CGMA_INIT_HDRS
  InitCGMA.hpp
)

add_library( cgma_init
  ${CGMA_INIT_SRCS}
)

target_link_libraries(cgma_init
  LINK_PUBLIC
    ${CGMA_LIBS}
  LINK_PRIVATE
    ${CGMA_DEPLIBS}
)
set(CGMA_LIBS "${CGMA_LIBS}" cgma_init PARENT_SCOPE)
set(CGMA_DEPLIBS "${CGMA_DEPLIBS}" PARENT_SCOPE)
add_definitions(${CGMA_ACIS_DEFINES} ${CGMA_OCC_DEFINES})

install(
  TARGETS cgma_init
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
)
install(
  FILES ${CGMA_INIT_HDRS}
  DESTINATION include
  COMPONENT Development
)
