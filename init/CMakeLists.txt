project( cgma_init )

include_directories(
  ${cubit_geom_SOURCE_DIR}
  ${cubit_util_SOURCE_DIR}
  ${cubit_geom_SOURCE_DIR}/ACIS
  ${cubit_geom_SOURCE_DIR}/virtual
  ${cubit_geom_SOURCE_DIR}/facet
  ${cubit_geom_SOURCE_DIR}/cholla
  ${cgma_init_SOURCE_DIR}
)

add_library( cgma_init
  InitCGMA.cpp
)

target_link_libraries( cgma_init
  cubit_geom
  cubit_util
)
