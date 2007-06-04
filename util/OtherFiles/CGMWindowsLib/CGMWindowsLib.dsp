# Microsoft Developer Studio Project File - Name="CGMWindowsLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=CGMWindowsLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "CGMWindowsLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "CGMWindowsLib.mak" CFG="CGMWindowsLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "CGMWindowsLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "CGMWindowsLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

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
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "../../../geom/Cholla" /I "../../../util" /I "../../../geom" /I "../../../geom/virtual" /I "../../../geom/facet" /D "NT" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D "MMGR_FREELIST" /D CUBIT_ACIS_VERSION=700 /D "ACIS_DLL" /FD /GZ /c
# SUBTRACT CPP /YX
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

# Name "CGMWindowsLib - Win32 Release"
# Name "CGMWindowsLib - Win32 Debug"
# Begin Group "geom"

# PROP Default_Filter ""
# Begin Group "ACIS"

# PROP Default_Filter ""
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisBridge.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisEdgeTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisFacetManager.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisHealerTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisModifyEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisQueryEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisSurfaceTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisTweakTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_cubit_owner.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_gtc.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_gtc_name.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_snl.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_snl_simple.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\BodyACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\CoEdgeACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\CurveACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\LoopACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\LumpACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\PointACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\ShellACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\SurfaceACIS.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/ACIS" /I "$(ACIS_DIR)\ag" /I "$(ACIS_DIR)\amain" /I "$(ACIS_DIR)\amfc" /I "$(ACIS_DIR)\blnd" /I "$(ACIS_DIR)\bool" /I "$(ACIS_DIR)\br" /I "$(ACIS_DIR)\clr" /I "$(ACIS_DIR)\covr" /I "$(ACIS_DIR)\cstr" /I "$(ACIS_DIR)\eulr" /I "$(ACIS_DIR)\fct" /I "$(ACIS_DIR)\ga" /I "$(ACIS_DIR)\gi" /I "$(ACIS_DIR)\gl" /I "$(ACIS_DIR)\heal" /I "$(ACIS_DIR)\iges" /I "$(ACIS_DIR)\ihl" /I "$(ACIS_DIR)\intr" /I "$(ACIS_DIR)\kern" /I "$(ACIS_DIR)\lop" /I "$(ACIS_DIR)\lopt" /I "$(ACIS_DIR)\base" /I "$(ACIS_DIR)\ofst" /I "$(ACIS_DIR)\oper" /I "$(ACIS_DIR)\part" /I "$(ACIS_DIR)\rbase" /I "$(ACIS_DIR)\rbi" /I "$(ACIS_DIR)\rem" /I "$(ACIS_DIR)\rom" /I "$(ACIS_DIR)\skin" /I "$(ACIS_DIR)\step" /I "$(ACIS_DIR)\swp" /I "$(ACIS_DIR)\trans" /I "$(ACIS_DIR)\law" /I "$(ACIS_DIR)\shl" /I "$(ACIS_DIR)\sbool" /I "$(ACIS_DIR)"

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisBridge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisEdgeTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisFacetManager.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisHealerTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisSurfaceTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\AcisTweakTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_cubit_owner.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_gtc.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_gtc_name.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_snl.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\attrib_snl_simple.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\BodyACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\CoEdgeACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\CurveACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\LoopACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\LumpACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\PointACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\ShellACIS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ACIS\SurfaceACIS.hpp
# End Source File
# End Group
# End Group
# Begin Group "virtual"

# PROP Default_Filter ""
# Begin Group "Source Files No. 1"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\virtual\CACompositeVG.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CAPartitionVG.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CAVirtualVG.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeBody.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeCurve.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeLump.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeSurface.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\Faceter.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\FaceterFacetData.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\FaceterPointData.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintBoundaryTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintPointData.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurve.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveGeometry.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveSegmented.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveVPSegGeo.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasitePoint.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionCurve.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionLump.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionManager.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionSurface.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartSurfTess.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PST_Data.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\SurfaceImprintTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGCoEdge.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGCurve.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGHiddenEntityManager.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGLoop.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGShell.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGSurface.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualGeometryEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualImprintTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualOSME.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VPFacetSurfGeom.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VPRelCoordSys.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files No. 1"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\virtual\CACompositeVG.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CAPartitionVG.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CAVirtualVG.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeBody.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeLump.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeSurface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\CompositeTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\Faceter.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\FaceterFacetData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\FaceterPointData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintBoundaryTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintLineSegment.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintMatchData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ImprintPointData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveGeometry.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveSegmented.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteCurveVPSegGeo.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasiteEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\ParasitePoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionLump.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionManager.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionSurface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartitionTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PartSurfTess.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\PST_Data.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\SurfaceImprintTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\TDVPRelCS.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGCoEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGHiddenEntityManager.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGLoop.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGShell.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VGSurface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualGeometryEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualImprintTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VirtualOSME.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VPFacetSurfGeom.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\virtual\VPRelCoordSys.hpp
# End Source File
# End Group
# End Group
# Begin Group "facet"

# PROP Default_Filter ""
# Begin Group "Source Files No. 2"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetBody.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetCoEdge.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetCurve.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetLoop.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetLump.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetModifyEngine.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetParamTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetPoint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetQueryEngine.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetShell.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetSurface.cpp
# End Source File
# End Group
# Begin Group "Header Files No. 2"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetBody.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetCoEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetCurveMesh.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetGeometryCreator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetLoop.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetLump.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetParamTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetPoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetPointMesh.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetShell.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facet\FacetSurface.hpp
# End Source File
# End Group
# End Group
# Begin Group "facetbool"

# PROP Default_Filter ""
# Begin Group "Source Files No. 3"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\facetbool\FacetedBoolean.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /I "../../../geom/facetbool"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facetbool\FSInt_stub.cpp
# End Source File
# End Group
# Begin Group "Header Files No. 3"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\facetbool\FacetedBoolean.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facetbool\FacetedSurfaceIntersector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\facetbool\FBDefines.hpp
# End Source File
# End Group
# End Group
# Begin Group "Source Files No. 4"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\AllocMemManagersGeom.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\AnalyticGeometryTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\AssemblySM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BasicTopologyEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Body.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BodySM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BoundingBoxTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BridgeManager.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAActuateSet.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CABodies.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CADeferredAttrib.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityColor.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityId.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityName.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAGroup.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAMergePartner.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAMergeStatus.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAUniqueId.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CGMApp.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Chain.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoEdge.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoEdgeSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoFace.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CollectionEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoVertex.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoVolume.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Csys.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CsysEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CsysOE.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitAttrib.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitAttribUser.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitPolygon.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitSimpleAttrib.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Curve.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CurveSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\DAG.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\DagDrawingTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\FacetTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\FastqClasses.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GDrawingTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryModifyEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryModifyTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryQueryEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryQueryTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryUtil.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeomMeasureTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GroupingEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GSaveOpen.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Loop.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\LoopSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Lump.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\LumpSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MergeTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MergeToolAssistant.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MidPlaneTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ModelEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ModelQueryEngine.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherRefEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherSolidModelEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ParamTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PartSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PlanarParamTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Point.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PointGridSearch.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PointSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefAssembly.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefCollection.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefCoordSys.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEdge.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntityFactory.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntityName.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefFace.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefGroup.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefPart.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefVertex.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefVolume.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SenseEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Shell.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ShellSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SplitSurfaceTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Surface.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfaceOverlapTool.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfaceSM.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfParamTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDCAGE.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDParallel.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDUniqueId.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TopologyBridge.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TopologyEntity.cpp

!IF  "$(CFG)" == "CGMWindowsLib - Win32 Release"

!ELSEIF  "$(CFG)" == "CGMWindowsLib - Win32 Debug"

# ADD CPP /GR

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files No. 4"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\AnalyticGeometryTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\AssemblySM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BasicTopologyEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Body.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BodySM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BoundingBoxTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\BridgeManager.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAActuateSet.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CABodies.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CADeferredAttrib.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityColor.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityId.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAEntityName.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAGroup.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAMergePartner.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAMergeStatus.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CAUniqueId.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CGMApp.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Chain.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoEdgeSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoFace.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CollectionEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoVertex.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CoVolume.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Csys.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CsysEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CsysOE.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitAttrib.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitAttribUser.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitPolygon.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CubitSimpleAttrib.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Curve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\CurveSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\DAG.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\DagDrawingTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\FacetEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\FacetTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\FastqClasses.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GDrawingTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryModifyTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryQueryTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeometryUtil.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GeomMeasureTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GroupingEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\GSaveOpen.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Loop.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\LoopSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Lump.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\LumpSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MergeTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MergeToolAssistant.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\MidPlaneTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ModelEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ModelQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherRefEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\OtherSolidModelEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ParamTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PartSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PlanarParamTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Point.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PointGridSearch.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\PointSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefAssembly.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefCollection.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefCoordSys.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntityFactory.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntityName.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefEntityNameMap.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefFace.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefGroup.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefPart.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefVertex.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\RefVolume.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SenseEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Shell.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\ShellSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SplitSurfaceTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Surface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfaceOverlapTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfaceSM.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\SurfParamTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDCAGE.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDCompare.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDParallel.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TDUniqueId.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TopologyBridge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\TopologyEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\UnMergeEvent.hpp
# End Source File
# End Group
# Begin Group "Cholla"

# PROP Default_Filter ""
# Begin Group "Source Files No. 6"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\Cholla\AllocMemManagersCholla.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CDrawingTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\Cholla.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaCurve.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaEngine.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaSkinTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaSurface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetEdge.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetEdgeData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitPoint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitPointData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitQuadFacet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitQuadFacetData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CurveFacetEvalTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\debug.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\FacetEntity.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\FacetEvalTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDFacetBoundaryEdge.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDFacetBoundaryPoint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDGeomFacet.cpp
# End Source File
# End Group
# Begin Group "Header Files No. 6"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CDrawingTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\Cholla.h
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaCurve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaEngine.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaPoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaSkinTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\ChollaSurface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacet.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitFacetEdgeData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitPoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitPointData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitQuadFacet.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CubitQuadFacetData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\CurveFacetEvalTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\debug.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\FacetEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\FacetEvalTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDFacetBoundaryEdge.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDFacetBoundaryPoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\geom\Cholla\TDGeomFacet.hpp
# End Source File
# End Group
# End Group
# End Group
# Begin Group "util"

# PROP Default_Filter ""
# Begin Group "Source Files No. 5"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\AllocMemManagersList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\AppUtil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\ArrayBasedContainer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CpuTimer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cubit2DPoint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitBox.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitCollection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitColorTable.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitContainer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitEntity.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitMessage.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitObservable.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitObserver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitPlane.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitStack.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitString.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitTransformMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitUtil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitVector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\DIntArray.cpp
# End Source File
# Begin Source File

SOURCE=..\..\DLList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\DynamicArray.cpp
# End Source File
# Begin Source File

SOURCE=..\..\GetLongOpt.cpp
# End Source File
# Begin Source File

SOURCE=..\..\GMem.cpp
# End Source File
# Begin Source File

SOURCE=..\..\IntersectionTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\MemoryBlock.cpp
# End Source File
# Begin Source File

SOURCE=..\..\MemoryManager.cpp
# End Source File
# Begin Source File

SOURCE=..\..\ParamCubitPlane.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Queue.cpp
# End Source File
# Begin Source File

SOURCE=..\..\RandomMersenne.cpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\SLList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\StubProgressTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\TextProgressTool.cpp
# End Source File
# Begin Source File

SOURCE=..\..\ToolData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\ToolDataUser.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Tree.cpp
# End Source File
# End Group
# Begin Group "Header Files No. 5"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\AppUtil.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ArrayBasedContainer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CastTo.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CpuTimer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\Cubit2DPoint.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitBox.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitBoxStruct.h
# End Source File
# Begin Source File

SOURCE=..\..\CubitCollection.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitColorConstants.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitColorTable.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitContainer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitDefines.h
# End Source File
# Begin Source File

SOURCE=..\..\CubitEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitEvent.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitInputFile.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitObservable.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitObserver.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitPlane.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitPlaneStruct.h
# End Source File
# Begin Source File

SOURCE=..\..\CubitStack.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitString.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitTransformMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitUtil.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\CubitVectorStruct.h
# End Source File
# Begin Source File

SOURCE=..\..\database.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DIntArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DLCubitEdgeList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DLCubitFaceList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DLIList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DLList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DLMRefEdgeList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DMRefFaceArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DoubleListItem.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DrawingToolDefines.h
# End Source File
# Begin Source File

SOURCE=..\..\DRefEdgeArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DRefFaceArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DRefVertexArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DynamicArray.hpp
# End Source File
# Begin Source File

SOURCE=..\..\DynamicDLIIterator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ElementType.h
# End Source File
# Begin Source File

SOURCE=..\..\FixedColorDefinitions.h
# End Source File
# Begin Source File

SOURCE=..\..\GeometryDefines.h
# End Source File
# Begin Source File

SOURCE=..\..\GetLongOpt.hpp
# End Source File
# Begin Source File

SOURCE=..\..\GMem.hpp
# End Source File
# Begin Source File

SOURCE=..\..\GUIProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\IndexedDouble.hpp
# End Source File
# Begin Source File

SOURCE=..\..\IntersectionTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\InvalidEntity.hpp
# End Source File
# Begin Source File

SOURCE=..\..\LocalStart.h
# End Source File
# Begin Source File

SOURCE=..\..\MemoryBlock.hpp
# End Source File
# Begin Source File

SOURCE=..\..\MemoryManager.hpp
# End Source File
# Begin Source File

SOURCE=..\..\MergeEvent.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ParamCubitPlane.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\Queue.hpp
# End Source File
# Begin Source File

SOURCE=..\..\RTree.hpp
# End Source File
# Begin Source File

SOURCE=..\..\RTreeNode.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLCADeferredAttribList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLCAMergePartnerList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLCubitAttribList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLDoubleList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLHexList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLIndexedDoubleList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLMRefEdgeLengthList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLTDAutoDetailList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SDLUniqueIdList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\SLList.hpp
# End Source File
# Begin Source File

SOURCE=..\..\StubProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\TDCellIndex.hpp
# End Source File
# Begin Source File

SOURCE=..\..\TextProgressTool.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ToolData.hpp
# End Source File
# Begin Source File

SOURCE=..\..\ToolDataUser.hpp
# End Source File
# Begin Source File

SOURCE=..\..\Tree.hpp
# End Source File
# End Group
# End Group
# End Target
# End Project
