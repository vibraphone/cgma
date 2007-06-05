# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Override default defines with the ones we want from the configure script
DEFS = $(TEMPLATE_DEFS_INCLUDED) $(LITTLE_ENDIAN)

# The name of the library to build
lib_LTLIBRARIES = 
if BUILD_CGM
  lib_LTLIBRARIES += libcubit_util.la
endif

# The directory where headers will be installed
libcubit_util_la_includedir = $(includedir)

# The non-template sources
# If any of these do not need to be installed, move them
# to the _SOURCES list above.
libcubit_util_la_SOURCES = \
  AllocMemManagersList.cpp \
  AppUtil.cpp \
  ArrayBasedContainer.cpp \
  CpuTimer.cpp \
  Cubit2DPoint.cpp \
  CubitBox.cpp \
  CubitCollection.cpp \
  CubitContainer.cpp \
  CubitCoordinateSystem.cpp \
  CubitDynamicLoader.cpp \
  CubitEntity.cpp \
  CubitFile.cpp \
  CubitFileFEModel.cpp \
  CubitFileIOWrapper.cpp \
  CubitFileMetaData.cpp \
  CubitFileUtil.cpp \
  CubitMatrix.cpp \
  CubitMessage.cpp \
  CubitObservable.cpp \
  CubitObserver.cpp \
  CubitPlane.cpp \
  CubitStack.cpp \
  CubitString.cpp \
  CubitTransformMatrix.cpp \
  CubitUtil.cpp \
  CubitVector.cpp \
  DIntArray.cpp \
  DLList.cpp \
  DynamicArray.cpp \
  GMem.cpp \
  GetLongOpt.cpp \
  GfxDebug.cpp \
  IntersectionTool.cpp \
  MemoryBlock.cpp \
  MemoryManager.cpp \
  ParamCubitPlane.cpp \
  PlanarParamTool.cpp \
  Queue.cpp \
  RandomMersenne.cpp \
  SDLList.cpp \
  SettingHandler.cpp \
  SettingHolder.cpp \
  StubProgressTool.cpp \
  TDUPtr.cpp \
  TextProgressTool.cpp \
  ToolData.cpp \
  ToolDataUser.cpp \
  Tree.cpp \
  TtyProgressTool.cpp 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_util_la_include_HEADERS = \
  AbstractTree.hpp \
  AppUtil.hpp \
  ArrayBasedContainer.hpp \
  CastTo.hpp \
  CpuTimer.hpp \
  Cubit2DPoint.hpp \
  CubitBox.hpp \
  CubitBoxStruct.h \
  CubitCollection.hpp \
  CubitColorConstants.hpp \
  CubitContainer.hpp \
  CubitCoordinateSystem.hpp \
  CubitDefines.h \
  CubitDynamicLoader.hpp \
  CubitEntity.hpp \
  CubitEvent.hpp \
  CubitEventDefines.h \
  CubitFile.hpp \
  CubitFileFEModel.hpp \
  CubitFileIOWrapper.hpp \
  CubitFileMetaData.hpp \
  CubitFileUtil.hpp \
  CubitInputFile.hpp \
  CubitMatrix.hpp \
  CubitMessage.hpp \
  CubitObservable.hpp \
  CubitObserver.hpp \
  CubitPlane.hpp \
  CubitPlaneStruct.h \
  CubitStack.hpp \
  CubitString.hpp \
  CubitTransformMatrix.hpp \
  CubitUtil.hpp \
  CubitUtilConfigure.h \
  CubitVector.hpp \
  CubitVectorStruct.h \
  DIntArray.hpp \
  DLCubitEdgeList.hpp \
  DLCubitFaceList.hpp \
  DLIList.hpp \
  DLList.hpp \
  DLMRefEdgeList.hpp \
  DMRefFaceArray.hpp \
  DRefEdgeArray.hpp \
  DRefFaceArray.hpp \
  DRefVertexArray.hpp \
  DoubleListItem.hpp \
  DrawingToolDefines.h \
  DynamicArray.hpp \
  DynamicDLIIterator.hpp \
  DynamicTreeIterator.hpp \
  ElementType.h \
  GMem.hpp \
  GeometryDefines.h \
  GetLongOpt.hpp \
  GfxDebug.hpp \
  IGUIObservers.hpp \
  IdSetEvent.hpp \
  IndexedDouble.hpp \
  IntersectionTool.hpp \
  InvalidEntity.hpp \
  KDDTree.hpp \
  KDDTreeNode.hpp \
  LocalStart.h \
  MemoryBlock.hpp \
  MemoryManager.hpp \
  MergeEvent.hpp \
  OctTree.hpp \
  OctTreeCell.hpp \
  ParamCubitPlane.hpp \
  ParamTool.hpp \
  PlanarParamTool.hpp \
  PriorityQueue.hpp \
  ProgressTool.hpp \
  Queue.hpp \
  RStarTree.hpp \
  RStarTreeNode.hpp \
  RTree.hpp \
  RTreeNode.hpp \
  RandomMersenne.hpp \
  SDLCADeferredAttribList.hpp \
  SDLCAMergePartnerList.hpp \
  SDLCubitAttribList.hpp \
  SDLDoubleList.hpp \
  SDLHexList.hpp \
  SDLIndexedDoubleList.hpp \
  SDLList.hpp \
  SDLMRefEdgeLengthList.hpp \
  SDLTDAutoDetailList.hpp \
  SettingHandler.hpp \
  SettingHolder.hpp \
  StubProgressTool.hpp \
  TDCellIndex.hpp \
  TDExport.hpp \
  TDUPtr.hpp \
  TDVector.hpp \
  TextProgressTool.hpp \
  ToolData.hpp \
  ToolDataUser.hpp \
  Tree.hpp \
  TtyProgressTool.hpp \
  database.hpp

# If template defs are included, then the template definitions
# need to be installed with the headers.  Otherwise they need
# to be compiled.
if INCLUDE_TEMPLATE_DEFS
  libcubit_util_la_include_HEADERS += \
      DLIList.cpp \
      KDDTree.cpp \
      KDDTreeNode.cpp \
      OctTreeCell.cpp \
      OctTree.cpp \
      PriorityQueue.cpp \
      RTree.cpp \
      RTreeNode.cpp \
      RStarTree.cpp \
      RStarTreeNode.cpp
else
  libcubit_util_la_SOURCES += \
      DLIList.cpp \
      KDDTree.cpp \
      KDDTreeNode.cpp \
      OctTreeCell.cpp \
      OctTree.cpp \
      PriorityQueue.cpp \
      RTree.cpp \
      RTreeNode.cpp \
      RStarTree.cpp \
      RStarTreeNode.cpp
endif

