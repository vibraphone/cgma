# Microsoft Developer Studio Project File - Name="libcubit_facet" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libcubit_facet - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_facet.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_facet.mak" CFG="libcubit_facet - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libcubit_facet - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "libcubit_facet - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "libcubit_facet - Win32 Debug"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /Zi /Od /I "../../util" /I ".." /I "../Cholla" /I "../facetbool" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /FD /GZ /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libcubit_facet - Win32 Release"

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
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W2 /GR /GX /O2 /I "../../util" /I ".." /I "../Cholla" /I "../facetbool" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /FD /c
# SUBTRACT CPP /YX
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

# Name "libcubit_facet - Win32 Debug"
# Name "libcubit_facet - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter "c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\FacetAttrib.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetAttribSet.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetBody.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetboolInterface.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetCoEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetLoop.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetLump.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetModifyEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetParamTool.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetQueryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetShell.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\GridSearchTree.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\FacetAttrib.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetAttribSet.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetBody.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetCoEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetLoop.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetLump.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetParamTool.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetShell.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\GridSearchTree.hpp
# End Source File
# Begin Source File

SOURCE=.\GridSearchTreeNode.hpp
# End Source File
# End Group
# Begin Source File

SOURCE=.\Readme.txt
# End Source File
# End Target
# End Project
