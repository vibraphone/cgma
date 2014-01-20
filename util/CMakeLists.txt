project(cubit_util)

# User option to build with concurrent execution support.
option(BUILD_WITH_CONCURRENT_SUPPORT "Build Cubit with concurrent execution support" OFF)

set(CGMA_UTIL_SRCS
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
  CubitEventDispatcher.cpp
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

set(CGMA_UTIL_HDRS
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
  CubitEvent.hpp
  CubitEventDefines.h
  CubitEventDispatcher.hpp
  CubitFile.hpp
  CubitFileFEModel.hpp
  CubitFileIOWrapper.hpp
  CubitFileMetaData.hpp
  CubitFileSimModel.hpp
  CubitFileUtil.hpp
  CubitInputFile.hpp
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
  ${cubit_util_BINARY_DIR}/CubitUtilConfigure.h

  AbstractTree.hpp
  CastTo.hpp
  CubitBoxStruct.h
  CubitColorConstants.hpp
  CubitPlaneStruct.h
  CubitVectorStruct.h
  DMRefFaceArray.hpp
  DRefEdgeArray.hpp
  DRefFaceArray.hpp
  DRefVertexArray.hpp
  DoubleListItem.hpp
  DrawingToolDefines.h
  DynamicDLIIterator.hpp
  DynamicTreeIterator.hpp
  ElementType.h
  GeometryDefines.h
  IGUIObservers.hpp
  IdSetEvent.hpp
  IndexedDouble.hpp
  InvalidEntity.hpp
  KDDTree.cpp
  KDDTree.hpp
  KDDTreeNode.cpp
  KDDTreeNode.hpp
  LocalStart.h
  ManagedPtrVector.hpp
  MergeEvent.hpp
  OctTree.cpp
  OctTree.hpp
  OctTreeCell.cpp
  OctTreeCell.hpp
  ParamTool.hpp
  PriorityQueue.cpp
  PriorityQueue.hpp
  ProgressTool.hpp
  RStarTree.cpp
  RStarTree.hpp
  RStarTreeNode.cpp
  RStarTreeNode.hpp
  RTree.cpp
  RTree.hpp
  RTreeNode.cpp
  RTreeNode.hpp
  SDLCADeferredAttribList.hpp
  SDLCAMergePartnerList.hpp
  SDLCubitAttribList.hpp
  SDLDoubleList.hpp
  SDLHexList.hpp
  SDLIndexedDoubleList.hpp
  SDLMRefEdgeLengthList.hpp
  SDLTDAutoDetailList.hpp
  TDCellIndex.hpp
  TDVector.hpp
  database.hpp
)

if(BUILD_WITH_CONCURRENT_SUPPORT)
  set(CGMA_UTIL_SRCS
    ${CGMA_UTIL_SRCS}
    CubitConcurrentApi.cpp
    CubitQtConcurrentApi.cpp
  )
  set(CGMA_UTIL_HDRS
    ${CGMA_UTIL_HDRS}
    CubitConcurrentApi.h
    CubitQtConcurrentApi.h
  )
  find_package(Qt4 REQUIRED QtCore)
  include(${QT_USE_FILE})
endif()


include_directories(${cubit_util_BINARY_DIR})

set(CUBIT_UTIL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  ${cubit_util_SOURCE_DIR}/CubitUtilConfigure.h.in
  ${cubit_util_BINARY_DIR}/CubitUtilConfigure.h
  @ONLY)

add_library(cubit_util ${CGMA_UTIL_SRCS})
target_link_libraries(cubit_util ${CMAKE_DL_LIBS})

if(CUBIT_UTIL_NAME)
  set_target_properties(cubit_util
    PROPERTIES OUTPUT_NAME "${CUBIT_UTIL_NAME}")
endif()

if(WIN32)
  target_link_libraries(cubit_util shlwapi)
endif(WIN32)

if(CUBIT_LIBRARY_PROPERTIES)
  set_target_properties(cubit_util
    PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
endif(CUBIT_LIBRARY_PROPERTIES)

if(BUILD_WITH_CONCURRENT_SUPPORT)
  target_link_libraries(cubit_util ${QT_LIBRARIES})
endif(BUILD_WITH_CONCURRENT_SUPPORT)

if(NOT CMAKE_INSTALL_BINARY_DIR)
  set(CMAKE_INSTALL_BINARY_DIR bin)
endif(NOT CMAKE_INSTALL_BINARY_DIR)

if(NOT CMAKE_INSTALL_INCLUDE_DIR)
  set(CMAKE_INSTALL_INCLUDE_DIR include)
endif(NOT CMAKE_INSTALL_INCLUDE_DIR)

if(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)
  if(BUILD_SHARED_LIBS)
    install(
      TARGETS cubit_util ${CUBIT_UTIL_EXPORT_GROUP}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
    )
  endif(BUILD_SHARED_LIBS)
else(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)
  install(
    TARGETS cubit_util ${CUBIT_UTIL_EXPORT_GROUP}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Runtime
    ARCHIVE DESTINATION ${CMAKE_INSTALL_BINARY_DIR} COMPONENT Development
  )
  if(NOT CUBIT_UTIL_INSTALL_NO_DEVELOPMENT_HEADERS)
    install(
      FILES ${CGMA_UTIL_HDRS}
      DESTINATION ${CMAKE_INSTALL_INCLUDE_DIR}
      COMPONENT Development
    )
  endif(NOT CUBIT_UTIL_INSTALL_NO_DEVELOPMENT_HEADERS)
endif(CUBIT_UTIL_INSTALL_NO_DEVELOPMENT)
