# Microsoft Developer Studio Project File - Name="libcubit_ACIS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libcubit_ACIS - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_ACIS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_ACIS.mak" CFG="libcubit_ACIS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libcubit_ACIS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "libcubit_ACIS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libcubit_ACIS - Win32 Debug"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE F90 /include:"Debug/"
# ADD F90 /include:"Debug/"
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /Gf /I "../../../util" /I "../../../list" /I "../../../geom" /I "../../../" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "OTPRO" /D "CUBIT_GUI" /D "NT" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D "MMGR_FREELIST" /D CUBIT_ACIS_VERSION=600 /FD /GZ /c
# SUBTRACT BASE CPP /YX /Yc /Yu
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "../../util" /I ".." /I "../Cholla" /I "$(ACIS_DIR)/include" /D "_DEBUG" /D "NO_FLEXLM" /D "NO_USAGE_TRACKING" /D "NT_DLLD" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D CUBIT_ACIS_VERSION=1302 /D "_SECDLL" /D "ACIS_DLL" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libcubit_ACIS - Win32 Release"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE F90 /include:"Release/"
# ADD F90 /include:"Release/"
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W2 /GR /GX /O2 /I "../../util" /I ".." /I "../Cholla" /I "$(ACIS_DIR)/include" /D "NDEBUG" /D "NT_DLL" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D CUBIT_ACIS_VERSION=1302 /D "_SECDLL" /D "ACIS_DLL" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /FD /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libcubit_ACIS - Win32 Debug"
# Name "libcubit_ACIS - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\AcisBridge.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisEdgeTool.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisFacetManager.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisHealerTool.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisModifyEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisQueryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisSurfaceTool.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisToolUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\AcisTweakTool.cpp
# End Source File
# Begin Source File

SOURCE=.\attrib_cubit_owner.cpp
# End Source File
# Begin Source File

SOURCE=.\attrib_gtc.cpp
# End Source File
# Begin Source File

SOURCE=.\attrib_gtc_name.cpp
# End Source File
# Begin Source File

SOURCE=.\attrib_snl.cpp
# End Source File
# Begin Source File

SOURCE=.\attrib_snl_simple.cpp
# End Source File
# Begin Source File

SOURCE=.\BodyACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\CoEdgeACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\CurveACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\LoopACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\LumpACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\PointACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\ShellACIS.cpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceACIS.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\AcisBridge.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisEdgeTool.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisFacetManager.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisHealerTool.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisSurfaceTool.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisToolUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisTweakTool.hpp
# End Source File
# Begin Source File

SOURCE=.\AcisTypes.h
# End Source File
# Begin Source File

SOURCE=.\attrib_cubit_owner.hpp
# End Source File
# Begin Source File

SOURCE=.\attrib_gtc.hpp
# End Source File
# Begin Source File

SOURCE=.\attrib_gtc_name.hpp
# End Source File
# Begin Source File

SOURCE=.\attrib_snl.hpp
# End Source File
# Begin Source File

SOURCE=.\attrib_snl_simple.hpp
# End Source File
# Begin Source File

SOURCE=.\BodyACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\CoEdgeACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\CurveACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\decl_none.h
# End Source File
# Begin Source File

SOURCE=.\LoopACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\LumpACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\PointACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\ShellACIS.hpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceACIS.hpp
# End Source File
# End Group
# End Target
# End Project
