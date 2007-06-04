# Microsoft Developer Studio Project File - Name="libcubit_virtual" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libcubit_virtual - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_virtual.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_virtual.mak" CFG="libcubit_virtual - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libcubit_virtual - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "libcubit_virtual - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libcubit_virtual - Win32 Debug"

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
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /I "../../../util" /I "../../../list" /I "../../../geom" /I "../../../" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "OTPRO" /D "CUBIT_GUI" /D "NT" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D "MMGR_FREELIST" /D CUBIT_ACIS_VERSION=503 /FD /GZ /c
# SUBTRACT BASE CPP /YX /Yc /Yu
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /Zi /Od /Gf /I "../facetbool" /I "../../util" /I ".." /I "../Cholla" /I "../facet" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_SECDLL" /D "NO_FLEXLM" /D "NO_USAGE_TRACKING" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /FD /GZ /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libcubit_virtual - Win32 Release"

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
# ADD CPP /nologo /MD /W2 /GR /GX /O2 /I "../../util" /I ".." /I "../Cholla" /I "../facet" /I "../facetbool" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_SECDLL" /D "NO_FLEXLM" /D "NO_USAGE_TRACKING" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /YX"stdafx.h" /FD /c
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

# Name "libcubit_virtual - Win32 Debug"
# Name "libcubit_virtual - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AllocMemManagersVirtual.cpp
# End Source File
# Begin Source File

SOURCE=.\CACompositeVG.cpp
# End Source File
# Begin Source File

SOURCE=.\CAPartitionVG.cpp
# End Source File
# Begin Source File

SOURCE=.\CAVirtualVG.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeAttrib.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeBody.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCoEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCoSurf.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeGeom.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeLoop.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeLump.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositePoint.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeShell.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\CompositeTool.cpp
# End Source File
# Begin Source File

SOURCE=.\CompSurfFacets.cpp
# End Source File
# Begin Source File

SOURCE=.\Faceter.cpp
# End Source File
# Begin Source File

SOURCE=.\FaceterFacetData.cpp
# End Source File
# Begin Source File

SOURCE=.\FaceterPointData.cpp
# End Source File
# Begin Source File

SOURCE=.\FacetProjectTool.cpp
# End Source File
# Begin Source File

SOURCE=.\HiddenEntitySet.cpp
# End Source File
# Begin Source File

SOURCE=.\ImprintBoundaryTool.cpp
# End Source File
# Begin Source File

SOURCE=.\ImprintLineSegment.cpp
# End Source File
# Begin Source File

SOURCE=.\ImprintPointData.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionBody.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCoEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCoSurf.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLoop.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLump.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLumpImprint.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionShell.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\PartitionTool.cpp
# End Source File
# Begin Source File

SOURCE=.\PartPTCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\PartSurfFacetTool.cpp
# End Source File
# Begin Source File

SOURCE=.\PST_Data.cpp
# End Source File
# Begin Source File

SOURCE=.\SegmentedCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\SubCurve.cpp
# End Source File
# Begin Source File

SOURCE=.\SubEntitySet.cpp
# End Source File
# Begin Source File

SOURCE=.\SubSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\TDVGFacetOwner.cpp
# End Source File
# Begin Source File

SOURCE=.\TDVGFacetSplit.cpp
# End Source File
# Begin Source File

SOURCE=.\VGArray.cpp
# End Source File
# Begin Source File

SOURCE=.\VGLoopTool.cpp
# End Source File
# Begin Source File

SOURCE=.\VirtualGeometryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\VirtualImprintTool.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\CACompositeVG.hpp
# End Source File
# Begin Source File

SOURCE=.\CAPartitionVG.hpp
# End Source File
# Begin Source File

SOURCE=.\CAVirtualVG.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeAttrib.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeBody.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCoEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCoSurf.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeGeom.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeLoop.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeLump.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositePoint.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeShell.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeTool.hpp
# End Source File
# Begin Source File

SOURCE=.\CompSurfFacets.hpp
# End Source File
# Begin Source File

SOURCE=.\Faceter.hpp
# End Source File
# Begin Source File

SOURCE=.\FaceterFacetData.hpp
# End Source File
# Begin Source File

SOURCE=.\FaceterPointData.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetProjectTool.hpp
# End Source File
# Begin Source File

SOURCE=.\HiddenEntitySet.hpp
# End Source File
# Begin Source File

SOURCE=.\ImprintBoundaryTool.hpp
# End Source File
# Begin Source File

SOURCE=.\ImprintLineSegment.hpp
# End Source File
# Begin Source File

SOURCE=.\ImprintMatchData.hpp
# End Source File
# Begin Source File

SOURCE=.\ImprintPointData.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionBody.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCoEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCoSurf.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLoop.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLump.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionLumpImprint.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionShell.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\PartitionTool.hpp
# End Source File
# Begin Source File

SOURCE=.\PartPTCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\PST_Data.hpp
# End Source File
# Begin Source File

SOURCE=.\SegmentedCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\SubCurve.hpp
# End Source File
# Begin Source File

SOURCE=.\SubEntitySet.hpp
# End Source File
# Begin Source File

SOURCE=.\SubSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\TDVGFacetOwner.hpp
# End Source File
# Begin Source File

SOURCE=.\TDVGFacetSplit.hpp
# End Source File
# Begin Source File

SOURCE=.\VGArray.hpp
# End Source File
# Begin Source File

SOURCE=.\VGDefines.h
# End Source File
# Begin Source File

SOURCE=.\VGLoopTool.hpp
# End Source File
# Begin Source File

SOURCE=.\VirtualGeometryEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\VirtualImprintTool.hpp
# End Source File
# End Group
# End Target
# End Project
