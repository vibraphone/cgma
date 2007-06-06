
PROJECT(gtcAttrib)

# Disallow in-source build
STRING(COMPARE EQUAL "${gtcAttrib_SOURCE_DIR}"
 "${gtcAttrib_BINARY_DIR}" INSOURCE)
IF(INSOURCE)
 MESSAGE(FATAL_ERROR "gtcAttrib requires an out of source Build. Please create a separate binary directory and run CMake there.")
 ENDIF(INSOURCE)

SET(GTC_ATTRIB_SRCS
    src/attrib_edge_parametric_length.C 
    src/attrib_gtc.C
    src/attrib_mark.C
    src/attrib_name.C
    )

# add acis includes to current project
INCLUDE_DIRECTORIES(${ACIS_INCLUDE_DIR} ${gtcAttrib_SOURCE_DIR}/incl)

#need to define ACIS_VER since it's used in the source here 
ADD_DEFINITIONS(-DACIS_VER=${ACIS_VERSION})

# build cubit acis library from sources
LINK_DIRECTORIES(${ACIS_LIBRARY_DIR} ${ACIS_LIBRARY_DIR_DEBUG})
ADD_LIBRARY(gtcAttrib SHARED ${GTC_ATTRIB_SRCS})

TARGET_LINK_LIBRARIES(gtcAttrib ${ACIS_LIBS})


IF(UNIX)
ADD_CUSTOM_TARGET(switch_acis ${CMAKE_MAKE_PROGRAM} clean)
ENDIF(UNIX)

INSTALL_TARGETS(/${CMAKE_INSTALL_BINARY_DIR} gtcAttrib)

