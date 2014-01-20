project(CGMA)
cmake_minimum_required(VERSION 2.8)

add_subdirectory(util)
add_subdirectory(init)
add_subdirectory(geom)
add_subdirectory(itaps)

install(
  FILES cgm_version.h
  DESTINATION include
  COMPONENT Development
)
