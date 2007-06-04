# Microsoft Developer Studio Project File - Name="cholla" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cholla - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cholla.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cholla.mak" CFG="cholla - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cholla - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cholla - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cholla - Win32 Release"

# PROP BASE Use_MFC 0
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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "../../util" /D "NDEBUG" /D "NT" /D "WIN32" /D "_MBCS" /D "_LIB" /D "TEMPLATE_DEFS_INCLUDED" /D "_AFXDLL" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "cholla - Win32 Debug"

# PROP BASE Use_MFC 0
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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "../../util" /D "_DEBUG" /D "NT" /D "WIN32" /D "_MBCS" /D "_LIB" /D "TEMPLATE_DEFS_INCLUDED" /D "_AFXDLL" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "cholla - Win32 Release"
# Name "cholla - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AllocMemManagersCholla.cpp
# End Source File
# Begin Source File

SOURCE=.\Cholla.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaSkinTool.cpp
# End Source File
# Begin Source File

SOURCE=.\ChollaSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacet.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetData.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetEdgeData.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitPointData.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitQuadFacet.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitQuadFacetData.cpp
# End Source File
# Begin Source File

SOURCE=.\CurveFacetEvalTool.cpp
# End Source File
# Begin Source File

SOURCE=.\debug.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetDataUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetEvalTool.cpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetboolData.cpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetBoundaryEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetBoundaryPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\TDGeomFacet.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\ChollaCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\ChollaEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\ChollaEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\ChollaPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\ChollaSkinTool.hpp
# End Source File
# Begin Source File

SOURCE=.\ChollaSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacet.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetData.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitFacetEdgeData.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitPointData.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitQuadFacet.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitQuadFacetData.hpp
# End Source File
# Begin Source File

SOURCE=.\CurveFacetEvalTool.hpp
# End Source File
# Begin Source File

SOURCE=.\debug.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetDataUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetEvalTool.hpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetboolData.hpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetBoundaryEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\TDFacetBoundaryPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\TDGeomFacet.hpp
# End Source File
# End Group
# End Target
# End Project
