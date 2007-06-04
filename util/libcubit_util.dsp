# Microsoft Developer Studio Project File - Name="libcubit_util" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libcubit_util - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_util.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_util.mak" CFG="libcubit_util - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libcubit_util - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libcubit_util - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libcubit_util - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "NT" /D "TEMPLATE_DEFS_INCLUDED" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libcubit_util - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "NT" /D "TEMPLATE_DEFS_INCLUDED" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libcubit_util - Win32 Release"
# Name "libcubit_util - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AllocMemManagersList.cpp
# End Source File
# Begin Source File

SOURCE=.\AppUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\ArrayBasedContainer.cpp
# End Source File
# Begin Source File

SOURCE=.\CpuTimer.cpp
# End Source File
# Begin Source File

SOURCE=.\Cubit2DPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitBox.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitCollection.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitContainer.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFile.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileFEModel.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileIOWrapper.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileMetaData.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitMatrix.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitMessage.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitObservable.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitObserver.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitPlane.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitStack.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitString.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitTransformMatrix.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitVector.cpp
# End Source File
# Begin Source File

SOURCE=.\DIntArray.cpp
# End Source File
# Begin Source File

SOURCE=.\DLList.cpp
# End Source File
# Begin Source File

SOURCE=.\DynamicArray.cpp
# End Source File
# Begin Source File

SOURCE=.\GetLongOpt.cpp
# End Source File
# Begin Source File

SOURCE=.\GfxDebug.cpp
# End Source File
# Begin Source File

SOURCE=.\GMem.cpp
# End Source File
# Begin Source File

SOURCE=.\IntersectionTool.cpp
# End Source File
# Begin Source File

SOURCE=.\MemoryBlock.cpp
# End Source File
# Begin Source File

SOURCE=.\MemoryManager.cpp
# End Source File
# Begin Source File

SOURCE=.\ParamCubitPlane.cpp
# End Source File
# Begin Source File

SOURCE=.\Queue.cpp
# End Source File
# Begin Source File

SOURCE=.\RandomMersenne.cpp
# End Source File
# Begin Source File

SOURCE=.\SDLList.cpp
# End Source File
# Begin Source File

SOURCE=.\SettingHandler.cpp
# End Source File
# Begin Source File

SOURCE=.\SettingHolder.cpp
# End Source File
# Begin Source File

SOURCE=.\SLList.cpp
# End Source File
# Begin Source File

SOURCE=.\StubProgressTool.cpp
# End Source File
# Begin Source File

SOURCE=.\TDUPtr.cpp
# End Source File
# Begin Source File

SOURCE=.\TextProgressTool.cpp
# End Source File
# Begin Source File

SOURCE=.\ToolData.cpp
# End Source File
# Begin Source File

SOURCE=.\ToolDataUser.cpp
# End Source File
# Begin Source File

SOURCE=.\Tree.cpp
# End Source File
# Begin Source File

SOURCE=.\TtyProgressTool.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\AbstractTree.hpp
# End Source File
# Begin Source File

SOURCE=.\AppUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\ArrayBasedContainer.hpp
# End Source File
# Begin Source File

SOURCE=.\CastTo.hpp
# End Source File
# Begin Source File

SOURCE=.\CpuTimer.hpp
# End Source File
# Begin Source File

SOURCE=.\Cubit2DPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitBox.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitBoxStruct.h
# End Source File
# Begin Source File

SOURCE=.\CubitCollection.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitColorConstants.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitColorTable.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitContainer.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitDefines.h
# End Source File
# Begin Source File

SOURCE=.\CubitEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitEvent.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitEventDefines.h
# End Source File
# Begin Source File

SOURCE=.\CubitFile.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileFEModel.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileIOWrapper.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFileMetaData.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitInputFile.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitMatrix.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitMessage.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitObservable.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitObserver.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitPlane.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitPlaneStruct.h
# End Source File
# Begin Source File

SOURCE=.\CubitStack.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitString.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitTransformMatrix.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitVector.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitVectorStruct.h
# End Source File
# Begin Source File

SOURCE=.\database.hpp
# End Source File
# Begin Source File

SOURCE=.\DIntArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DLCubitEdgeList.hpp
# End Source File
# Begin Source File

SOURCE=.\DLCubitFaceList.hpp
# End Source File
# Begin Source File

SOURCE=.\DLIList.hpp
# End Source File
# Begin Source File

SOURCE=.\DLList.hpp
# End Source File
# Begin Source File

SOURCE=.\DLMRefEdgeList.hpp
# End Source File
# Begin Source File

SOURCE=.\DMRefFaceArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DoubleListItem.hpp
# End Source File
# Begin Source File

SOURCE=.\DrawingToolDefines.h
# End Source File
# Begin Source File

SOURCE=.\DRefEdgeArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DRefFaceArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DRefVertexArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DynamicArray.hpp
# End Source File
# Begin Source File

SOURCE=.\DynamicDLIIterator.hpp
# End Source File
# Begin Source File

SOURCE=.\ElementType.h
# End Source File
# Begin Source File

SOURCE=.\FixedColorDefinitions.h
# End Source File
# Begin Source File

SOURCE=.\GeometryDefines.h
# End Source File
# Begin Source File

SOURCE=.\GetLongOpt.hpp
# End Source File
# Begin Source File

SOURCE=.\GfxDebug.hpp
# End Source File
# Begin Source File

SOURCE=.\GMem.hpp
# End Source File
# Begin Source File

SOURCE=.\GUIProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=.\IGUIObservers.hpp
# End Source File
# Begin Source File

SOURCE=.\IndexedDouble.hpp
# End Source File
# Begin Source File

SOURCE=.\IntersectionTool.hpp
# End Source File
# Begin Source File

SOURCE=.\InvalidEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\KDDTree.hpp
# End Source File
# Begin Source File

SOURCE=.\KDDTreeNode.hpp
# End Source File
# Begin Source File

SOURCE=.\LocalStart.h
# End Source File
# Begin Source File

SOURCE=.\MemoryBlock.hpp
# End Source File
# Begin Source File

SOURCE=.\MemoryManager.hpp
# End Source File
# Begin Source File

SOURCE=.\MergeEvent.hpp
# End Source File
# Begin Source File

SOURCE=.\OctTree.hpp
# End Source File
# Begin Source File

SOURCE=.\OctTreeCell.hpp
# End Source File
# Begin Source File

SOURCE=.\ParamCubitPlane.hpp
# End Source File
# Begin Source File

SOURCE=.\PriorityQueue.hpp
# End Source File
# Begin Source File

SOURCE=.\ProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=.\Queue.hpp
# End Source File
# Begin Source File

SOURCE=.\RandomMersenne.hpp
# End Source File
# Begin Source File

SOURCE=.\RStarTree.hpp
# End Source File
# Begin Source File

SOURCE=.\RStarTreeNode.hpp
# End Source File
# Begin Source File

SOURCE=.\RTree.hpp
# End Source File
# Begin Source File

SOURCE=.\RTreeNode.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLCADeferredAttribList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLCAMergePartnerList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLCubitAttribList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLDoubleList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLHexList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLIndexedDoubleList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLMRefEdgeLengthList.hpp
# End Source File
# Begin Source File

SOURCE=.\SDLTDAutoDetailList.hpp
# End Source File
# Begin Source File

SOURCE=.\SettingHandler.hpp
# End Source File
# Begin Source File

SOURCE=.\SettingHolder.hpp
# End Source File
# Begin Source File

SOURCE=.\SLList.hpp
# End Source File
# Begin Source File

SOURCE=.\StubProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=.\TDCellIndex.hpp
# End Source File
# Begin Source File

SOURCE=.\TDUPtr.hpp
# End Source File
# Begin Source File

SOURCE=.\TextProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=.\ToolData.hpp
# End Source File
# Begin Source File

SOURCE=.\ToolDataUser.hpp
# End Source File
# Begin Source File

SOURCE=.\Tree.hpp
# End Source File
# Begin Source File

SOURCE=.\TtyProgressTool.hpp
# End Source File
# End Group
# End Target
# End Project
