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
  CommandFeedback.hpp \
  CpuTimer.cpp \
  Cubit2DPoint.cpp \
  CubitBox.cpp \
  CubitCollection.cpp \
  CubitContainer.cpp \
  CubitCoordinateSystem.cpp \
  CubitDynamicLoader.cpp \
  CubitEntity.cpp \
  CubitEventDispatcher.cpp \
  CubitFile.cpp \
  CubitFileFEModel.cpp \
  CubitFileIOWrapper.cpp \
  CubitFileMetaData.cpp \
  CubitFileSimModel.cpp \
  CubitFileUtil.cpp \
  CubitMatrix.cpp \
  CubitMessage.cpp \
  CubitObservable.cpp \
  CubitObserver.cpp \
  CubitPlane.cpp \
  CubitSparseMatrix.cpp \
  CubitStack.cpp \
  CubitString.cpp \
  CubitTransformMatrix.cpp \
  CubitUndo.cpp \
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

#if WITH_CUBIT
  libcubit_util_la_SOURCES += GlobalCommandFeedback.cpp \
      CubitProcess.cpp \
      CubitInstrumentation.cpp
#endif 

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
  CubitEventDispatcher.hpp \
  CubitEventDefines.h \
  CubitFile.hpp \
  CubitFileFEModel.hpp \
  CubitFileIOWrapper.hpp \
  CubitFileMetaData.hpp \
  CubitFileSimModel.hpp \
  CubitLoops.hpp \
  CubitUndo.hpp \
  CubitFileUtil.hpp \
  CubitInputFile.hpp \
  CubitMatrix.hpp \
  CubitMessage.hpp \
  CubitObservable.hpp \
  CubitObserver.hpp \
  CubitPlane.hpp \
  CubitPlaneStruct.h \
  CubitSparseMatrix.hpp \
  CubitStack.hpp \
  CubitString.hpp \
  CubitTransformMatrix.hpp \
  CubitUtil.hpp \
  CubitUtilConfigure.h \
  CubitVector.hpp \
  CubitVectorStruct.h \
  DIntArray.hpp \
  DLIList.hpp \
  DLList.hpp \
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
  FacetShapeDefs.hpp \
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
  ManagedPtrVector.hpp \
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
  TDUPtr.hpp \
  TDVector.hpp \
  TextProgressTool.hpp \
  ToolData.hpp \
  ToolDataUser.hpp \
  Tree.hpp \
  TtyProgressTool.hpp \
  database.hpp

#if WITH_CUBIT
  libcubit_util_la_include_HEADERS += GlobalCommandFeedback.hpp \
           CubitProcess.hpp \
           CubitInstrumentation.hpp
#endif 
 
# If template defs are included, then the template definitions
# need to be installed with the headers.  Otherwise they need
# to be compiled.
if INCLUDE_TEMPLATE_DEFS
  libcubit_util_la_include_HEADERS += \
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

libcubit_util_la_LIBADD = -ldl 
