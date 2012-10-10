
PROJECT(cubit_util)

# User option to build with concurrent execution support.
OPTION(BUILD_WITH_CONCURRENT_SUPPORT "Build Cubit with concurrent execution support" OFF)

SET(UTIL_SRCS
    AllocMemManagersList.cpp
    AppUtil.cpp
    ArrayBasedContainer.cpp
    CpuTimer.cpp
    Cubit2DPoint.cpp
    CubitBox.cpp
    CubitCollection.cpp
    CubitContainer.cpp
    CubitCoordinateSystem.cpp
    CubitDynamicLoader.cpp
    CubitEntity.cpp
    CubitFile.cpp
    CubitFileFEModel.cpp
    CubitFileIOWrapper.cpp
    CubitFileMetaData.cpp
    CubitFileSimModel.cpp
    CubitFileUtil.cpp
    CubitInstrumentation.cpp
    CubitMatrix.cpp
    CubitMessage.cpp
    CubitObservable.cpp
    CubitObserver.cpp
    CubitPlane.cpp
    CubitProcess.cpp
    CubitSparseMatrix.cpp
    CubitStack.cpp
    CubitString.cpp
    CubitTransformMatrix.cpp
    CubitUndo.cpp
    CubitUtil.cpp
    CubitVector.cpp
    DIntArray.cpp
    DLList.cpp
    DynamicArray.cpp
    GetLongOpt.cpp
    GfxDebug.cpp
    GlobalCommandFeedback.cpp
    GMem.cpp
    IntersectionTool.cpp
    MemoryBlock.cpp
    MemoryManager.cpp
    ParamCubitPlane.cpp
    PlanarParamTool.cpp
    Queue.cpp
    RandomMersenne.cpp
    SDLList.cpp
    SettingHandler.cpp
    SettingHolder.cpp
    StubProgressTool.cpp
    TDUPtr.cpp
    TextProgressTool.cpp
    ToolData.cpp
    ToolDataUser.cpp
    Tree.cpp
    TtyProgressTool.cpp
    )

SET(UTIL_HDRS
    AppUtil.hpp
    ArrayBasedContainer.hpp
    CpuTimer.hpp
    CommandFeedback.hpp
    Cubit2DPoint.hpp
    CubitBox.hpp
    CubitCollection.hpp
    CubitContainer.hpp
    CubitCoordinateSystem.hpp
    CubitDefines.h
    CubitDynamicLoader.hpp
    CubitEntity.hpp
    CubitEquation.hpp
    CubitFile.hpp
    CubitFileFEModel.hpp
    CubitFileIOWrapper.hpp
    CubitFileMetaData.hpp
    CubitFileSimModel.hpp
    CubitFileUtil.hpp
    CubitInstrumentation.hpp
    CubitLoops.hpp
    CubitMatrix.hpp
    CubitMessage.hpp
    CubitObservable.hpp
    CubitObserver.hpp
    CubitPlane.hpp
    CubitProcess.hpp
    CubitSparseMatrix.hpp
    CubitStack.hpp
    CubitString.hpp
    CubitTransformMatrix.hpp
    CubitUndo.hpp
    CubitUtil.hpp
    CubitVector.hpp
    DIntArray.hpp
    DLIList.hpp
    DLList.hpp
    DynamicArray.hpp
    FacetShapeDefs.hpp
    GetLongOpt.hpp
    GfxDebug.hpp
    GlobalCommandFeedback.hpp
    GMem.hpp
    IdHasher.hpp
    IntersectionTool.hpp
    MemoryBlock.hpp
    MemoryManager.hpp
    OrderedSet.hpp
    OrderedMap.hpp
    ParamCubitPlane.hpp
    PlanarParamTool.hpp
    Queue.hpp
    RandomMersenne.hpp
    SDLList.hpp
    SettingHandler.hpp
    SettingHolder.hpp
    StubProgressTool.hpp
    TDUPtr.hpp
    TextProgressTool.hpp
    ToolData.hpp
    ToolDataUser.hpp
    Tree.hpp
    TtyProgressTool.hpp
    VariableArray.hpp
    CubitEventDefines.h
    ${cubit_util_BINARY_DIR}/CubitUtilConfigure.h
  )

IF(BUILD_WITH_CONCURRENT_SUPPORT)
  SET(UTIL_SRCS
      ${UTIL_SRCS}
      CubitConcurrentApi.cpp
      CubitQtConcurrentApi.cpp
      )
  SET(UTIL_HDRS
      ${UTIL_HDRS}
      CubitConcurrentApi.h
      CubitQtConcurrentApi.h
      )
  FIND_PACKAGE(Qt4 REQUIRED QtCore)
  INCLUDE(${QT_USE_FILE})
ENDIF(BUILD_WITH_CONCURRENT_SUPPORT)


INCLUDE_DIRECTORIES(${cubit_util_BINARY_DIR})

SET(CUBIT_UTIL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
CONFIGURE_FILE(${cubit_util_SOURCE_DIR}/CubitUtilConfigure.h.in
               ${cubit_util_BINARY_DIR}/CubitUtilConfigure.h
               @ONLY)

IF(NOT CUBIT_UTIL_NAME)
  SET(CUBIT_UTIL_NAME cubit_util)
ENDIF(NOT CUBIT_UTIL_NAME)

ADD_LIBRARY(${CUBIT_UTIL_NAME} ${UTIL_SRCS} ${UTIL_HDRS})
TARGET_LINK_LIBRARIES(${CUBIT_UTIL_NAME} ${CMAKE_DL_LIBS})

IF(CUBIT_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(${CUBIT_UTIL_NAME}
                        PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
ENDIF(CUBIT_LIBRARY_PROPERTIES)

IF(BUILD_WITH_CONCURRENT_SUPPORT)
  TARGET_LINK_LIBRARIES(${CUBIT_UTIL_NAME} ${QT_LIBRARIES})
ENDIF(BUILD_WITH_CONCURRENT_SUPPORT)

IF(NOT CMAKE_INSTALL_BINARY_DIR)
  SET(CMAKE_INSTALL_BINARY_DIR bin)
ENDIF(NOT CMAKE_INSTALL_BINARY_DIR)

IF(NOT CMAKE_INSTALL_INCLUDE_DIR)
  SET(CMAKE_INSTALL_INCLUDE_DIR include)
ENDIF(NOT CMAKE_INSTALL_INCLUDE_DIR)

IF(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)
  IF(BUILD_SHARED_LIBS)
    INSTALL(TARGETS ${CUBIT_UTIL_NAME} ${CUBIT_UTIL_EXPORT_GROUP}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
            LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime)
  ENDIF(BUILD_SHARED_LIBS)
ELSE(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)
INSTALL(TARGETS ${CUBIT_UTIL_NAME} ${CUBIT_UTIL_EXPORT_GROUP}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development)

IF(NOT CUBIT_UTIL_INSTALL_NO_DEVELOPMENT_HEADERS)
INSTALL(FILES ${UTIL_HDRS}
        DESTINATION ${CMAKE_INSTALL_INCLUDE_DIR} COMPONENT Development)
ENDIF(NOT CUBIT_UTIL_INSTALL_NO_DEVELOPMENT_HEADERS)

ENDIF(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)

