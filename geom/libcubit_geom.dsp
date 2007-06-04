# Microsoft Developer Studio Project File - Name="libcubit_geom" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libcubit_geom - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_geom.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libcubit_geom.mak" CFG="libcubit_geom - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libcubit_geom - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "libcubit_geom - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libcubit_geom - Win32 Debug"

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
# ADD F90 /browser /include:"Debug/"
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /I "../../../util" /I "../../../list" /I "../../../geom" /I "../../../" /I "../../../geom/virtual" /I "../../../geom/ACIS" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "OTPRO" /D "CUBIT_GUI" /D "NT" /D "ACIS_3D" /D "ACIS_HEALER" /D "ACIS_LOCAL_OPS" /D "ACIS_IGES_TRANSLATOR" /D "ACIS_STEP_TRANSLATOR" /D "MMGR_FREELIST" /D CUBIT_ACIS_VERSION=503 /YX"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /Zi /Od /Gf /I "../util" /I "./Cholla" /I "./facet" /I "./virtual" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_SECDLL" /D "NO_FLEXLM" /D "NO_USAGE_TRACKING" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /YX"stdafx.h" /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libcubit_geom - Win32 Release"

# PROP BASE Use_MFC 1
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
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W2 /GR /GX /O2 /I "../util" /I "./Cholla" /I "./facet" /I "./virtual" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "NT" /D "_SECDLL" /D "NO_FLEXLM" /D "NO_USAGE_TRACKING" /D "_AFXDLL" /D "TEMPLATE_DEFS_INCLUDED" /FR /YX"stdafx.h" /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libcubit_geom - Win32 Debug"
# Name "libcubit_geom - Win32 Release"
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\AnalyticGeometryTool.hpp
# End Source File
# Begin Source File

SOURCE=.\AssemblySM.hpp
# End Source File
# Begin Source File

SOURCE=.\BasicTopologyEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\Body.hpp
# End Source File
# Begin Source File

SOURCE=.\BodySM.hpp
# End Source File
# Begin Source File

SOURCE=.\BoundaryConstrainTool.hpp
# End Source File
# Begin Source File

SOURCE=.\BoundingBoxTool.hpp
# End Source File
# Begin Source File

SOURCE=.\BridgeManager.hpp
# End Source File
# Begin Source File

SOURCE=.\CAActuateSet.hpp
# End Source File
# Begin Source File

SOURCE=.\CABodies.hpp
# End Source File
# Begin Source File

SOURCE=.\CADeferredAttrib.hpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityColor.hpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityId.hpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityName.hpp
# End Source File
# Begin Source File

SOURCE=.\CAGroup.hpp
# End Source File
# Begin Source File

SOURCE=.\CAMergePartner.hpp
# End Source File
# Begin Source File

SOURCE=.\CAMergeStatus.hpp
# End Source File
# Begin Source File

SOURCE=.\CAUniqueId.hpp
# End Source File
# Begin Source File

SOURCE=.\CGMApp.hpp
# End Source File
# Begin Source File

SOURCE=.\Chain.hpp
# End Source File
# Begin Source File

SOURCE=.\ChordalAxis.hpp
# End Source File
# Begin Source File

SOURCE=.\CoEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\CoEdgeSM.hpp
# End Source File
# Begin Source File

SOURCE=.\CoFace.hpp
# End Source File
# Begin Source File

SOURCE=.\CollectionEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\CompositeCombineEvent.hpp
# End Source File
# Begin Source File

SOURCE=.\CoVertex.hpp
# End Source File
# Begin Source File

SOURCE=.\CoVolume.hpp
# End Source File
# Begin Source File

SOURCE=.\Csys.hpp
# End Source File
# Begin Source File

SOURCE=.\CsysEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\CsysOE.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitAttrib.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitAttribUser.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitEvaluator.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitPolygon.hpp
# End Source File
# Begin Source File

SOURCE=.\CubitSimpleAttrib.hpp
# End Source File
# Begin Source File

SOURCE=.\Curve.hpp
# End Source File
# Begin Source File

SOURCE=.\CurveSM.hpp
# End Source File
# Begin Source File

SOURCE=.\CylinderEvaluator.hpp
# End Source File
# Begin Source File

SOURCE=.\DAG.hpp
# End Source File
# Begin Source File

SOURCE=.\DagDrawingTool.hpp
# End Source File
# Begin Source File

SOURCE=.\DagType.hpp
# End Source File
# Begin Source File

SOURCE=.\FacetorTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GDrawingTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GeomDataObserver.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryHealerEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryHealerTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryModifyEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryModifyTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryQueryTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GeometryUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\GeomMeasureTool.hpp
# End Source File
# Begin Source File

SOURCE=.\GeomPoint.hpp
# End Source File
# Begin Source File

SOURCE=.\GeomSeg.hpp
# End Source File
# Begin Source File

SOURCE=.\GeoNode.hpp
# End Source File
# Begin Source File

SOURCE=.\GeoTet.hpp
# End Source File
# Begin Source File

SOURCE=.\GroupingEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\GSaveOpen.hpp
# End Source File
# Begin Source File

SOURCE=.\IntermediateGeomEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\Loop.hpp
# End Source File
# Begin Source File

SOURCE=.\LoopParamTool.hpp
# End Source File
# Begin Source File

SOURCE=.\LoopSM.hpp
# End Source File
# Begin Source File

SOURCE=.\Lump.hpp
# End Source File
# Begin Source File

SOURCE=.\LumpSM.hpp
# End Source File
# Begin Source File

SOURCE=.\MedialTool2D.hpp
# End Source File
# Begin Source File

SOURCE=.\MedialTool3D.hpp
# End Source File
# Begin Source File

SOURCE=.\MergeTool.hpp
# End Source File
# Begin Source File

SOURCE=.\MergeToolAssistant.hpp
# End Source File
# Begin Source File

SOURCE=.\MidPlaneTool.hpp
# End Source File
# Begin Source File

SOURCE=.\ModelEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\ModelQueryEngine.hpp
# End Source File
# Begin Source File

SOURCE=.\OtherEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\OtherRefEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\ParamTool.hpp
# End Source File
# Begin Source File

SOURCE=.\PartSM.hpp
# End Source File
# Begin Source File

SOURCE=.\PlanarParamTool.hpp
# End Source File
# Begin Source File

SOURCE=.\Point.hpp
# End Source File
# Begin Source File

SOURCE=.\PointGridSearch.hpp
# End Source File
# Begin Source File

SOURCE=.\PointLoopFacetor.hpp
# End Source File
# Begin Source File

SOURCE=.\PointSM.hpp
# End Source File
# Begin Source File

SOURCE=.\RefAssembly.hpp
# End Source File
# Begin Source File

SOURCE=.\RefCollection.hpp
# End Source File
# Begin Source File

SOURCE=.\RefCoordSys.hpp
# End Source File
# Begin Source File

SOURCE=.\RefEdge.hpp
# End Source File
# Begin Source File

SOURCE=.\RefEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\RefEntityFactory.hpp
# End Source File
# Begin Source File

SOURCE=.\RefEntityName.hpp
# End Source File
# Begin Source File

SOURCE=.\RefEntityNameMap.hpp
# End Source File
# Begin Source File

SOURCE=.\RefFace.hpp
# End Source File
# Begin Source File

SOURCE=.\RefGroup.hpp
# End Source File
# Begin Source File

SOURCE=.\RefPart.hpp
# End Source File
# Begin Source File

SOURCE=.\RefVertex.hpp
# End Source File
# Begin Source File

SOURCE=.\RefVolume.hpp
# End Source File
# Begin Source File

SOURCE=.\SenseEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\Shell.hpp
# End Source File
# Begin Source File

SOURCE=.\ShellSM.hpp
# End Source File
# Begin Source File

SOURCE=.\SphereEvaluator.hpp
# End Source File
# Begin Source File

SOURCE=.\SplitSurfaceTool.hpp
# End Source File
# Begin Source File

SOURCE=.\Surface.hpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceOverlapFacet.hpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceOverlapTool.hpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceSM.hpp
# End Source File
# Begin Source File

SOURCE=.\SurfParamTool.hpp
# End Source File
# Begin Source File

SOURCE=.\TBOwner.hpp
# End Source File
# Begin Source File

SOURCE=.\TBOwnerSet.hpp
# End Source File
# Begin Source File

SOURCE=.\TDCAGE.hpp
# End Source File
# Begin Source File

SOURCE=.\TDChordal.hpp
# End Source File
# Begin Source File

SOURCE=.\TDCompare.hpp
# End Source File
# Begin Source File

SOURCE=.\TDDelaunay.hpp
# End Source File
# Begin Source File

SOURCE=.\TDInterpNode.hpp
# End Source File
# Begin Source File

SOURCE=.\TDParallel.hpp
# End Source File
# Begin Source File

SOURCE=.\TDSplitSurface.hpp
# End Source File
# Begin Source File

SOURCE=.\TDSurfaceOverlap.hpp
# End Source File
# Begin Source File

SOURCE=.\TDUniqueId.hpp
# End Source File
# Begin Source File

SOURCE=.\TDVector.hpp
# End Source File
# Begin Source File

SOURCE=.\TetFacetorTool.hpp
# End Source File
# Begin Source File

SOURCE=.\TopologyBridge.hpp
# End Source File
# Begin Source File

SOURCE=.\TopologyEntity.hpp
# End Source File
# Begin Source File

SOURCE=.\UnMergeEvent.hpp
# End Source File
# End Group
# Begin Group "geom"

# PROP Default_Filter "*.cpp"
# Begin Source File

SOURCE=.\AllocMemManagersGeom.cpp
# End Source File
# Begin Source File

SOURCE=.\AnalyticGeometryTool.cpp
# End Source File
# Begin Source File

SOURCE=.\AssemblySM.cpp
# End Source File
# Begin Source File

SOURCE=.\BasicTopologyEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\Body.cpp
# End Source File
# Begin Source File

SOURCE=.\BodySM.cpp
# End Source File
# Begin Source File

SOURCE=.\BoundingBoxTool.cpp
# End Source File
# Begin Source File

SOURCE=.\BridgeManager.cpp
# End Source File
# Begin Source File

SOURCE=.\CAActuateSet.cpp
# End Source File
# Begin Source File

SOURCE=.\CABodies.cpp
# End Source File
# Begin Source File

SOURCE=.\CADeferredAttrib.cpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityColor.cpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityId.cpp
# End Source File
# Begin Source File

SOURCE=.\CAEntityName.cpp
# End Source File
# Begin Source File

SOURCE=.\CAGroup.cpp
# End Source File
# Begin Source File

SOURCE=.\CAMergePartner.cpp
# End Source File
# Begin Source File

SOURCE=.\CAMergeStatus.cpp
# End Source File
# Begin Source File

SOURCE=.\CAUniqueId.cpp
# End Source File
# Begin Source File

SOURCE=.\CGMApp.cpp
# End Source File
# Begin Source File

SOURCE=.\Chain.cpp
# End Source File
# Begin Source File

SOURCE=.\ChordalAxis.cpp
# End Source File
# Begin Source File

SOURCE=.\CoEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\CoEdgeSM.cpp
# End Source File
# Begin Source File

SOURCE=.\CoFace.cpp
# End Source File
# Begin Source File

SOURCE=.\CollectionEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\CoVertex.cpp
# End Source File
# Begin Source File

SOURCE=.\CoVolume.cpp
# End Source File
# Begin Source File

SOURCE=.\Csys.cpp
# End Source File
# Begin Source File

SOURCE=.\CsysEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\CsysOE.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitAttrib.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitAttribManager.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitAttribUser.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitPolygon.cpp
# End Source File
# Begin Source File

SOURCE=.\CubitSimpleAttrib.cpp
# End Source File
# Begin Source File

SOURCE=.\Curve.cpp
# End Source File
# Begin Source File

SOURCE=.\CurveSM.cpp
# End Source File
# Begin Source File

SOURCE=.\CylinderEvaluator.cpp
# End Source File
# Begin Source File

SOURCE=.\DAG.cpp
# End Source File
# Begin Source File

SOURCE=.\DagDrawingTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GDrawingTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GeomDataObserver.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryHealerEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryHealerTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryModifyEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryModifyTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryQueryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryQueryTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GeometryUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\GeomMeasureTool.cpp
# End Source File
# Begin Source File

SOURCE=.\GeoNode.cpp
# End Source File
# Begin Source File

SOURCE=.\GeoTet.cpp
# End Source File
# Begin Source File

SOURCE=.\GroupingEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\GSaveOpen.cpp
# End Source File
# Begin Source File

SOURCE=.\Loop.cpp
# End Source File
# Begin Source File

SOURCE=.\LoopParamTool.cpp
# End Source File
# Begin Source File

SOURCE=.\LoopSM.cpp
# End Source File
# Begin Source File

SOURCE=.\Lump.cpp
# End Source File
# Begin Source File

SOURCE=.\LumpSM.cpp
# End Source File
# Begin Source File

SOURCE=.\MedialTool2D.cpp
# End Source File
# Begin Source File

SOURCE=.\MedialTool3D.cpp
# End Source File
# Begin Source File

SOURCE=.\MergeTool.cpp
# End Source File
# Begin Source File

SOURCE=.\MergeToolAssistant.cpp
# End Source File
# Begin Source File

SOURCE=.\MidPlaneTool.cpp
# End Source File
# Begin Source File

SOURCE=.\ModelEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\ModelQueryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\OldUnmergeCode.cpp
# End Source File
# Begin Source File

SOURCE=.\OtherEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\OtherRefEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\PartSM.cpp
# End Source File
# Begin Source File

SOURCE=.\PlanarParamTool.cpp
# End Source File
# Begin Source File

SOURCE=.\Point.cpp
# End Source File
# Begin Source File

SOURCE=.\PointGridSearch.cpp
# End Source File
# Begin Source File

SOURCE=.\PointLoopFacetor.cpp
# End Source File
# Begin Source File

SOURCE=.\PointSM.cpp
# End Source File
# Begin Source File

SOURCE=.\RefAssembly.cpp
# End Source File
# Begin Source File

SOURCE=.\RefCollection.cpp
# End Source File
# Begin Source File

SOURCE=.\RefCoordSys.cpp
# End Source File
# Begin Source File

SOURCE=.\RefEdge.cpp
# End Source File
# Begin Source File

SOURCE=.\RefEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\RefEntityFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\RefEntityName.cpp
# End Source File
# Begin Source File

SOURCE=.\RefFace.cpp
# End Source File
# Begin Source File

SOURCE=.\RefGroup.cpp
# End Source File
# Begin Source File

SOURCE=.\RefPart.cpp
# End Source File
# Begin Source File

SOURCE=.\RefVertex.cpp
# End Source File
# Begin Source File

SOURCE=.\RefVolume.cpp
# End Source File
# Begin Source File

SOURCE=.\SenseEntity.cpp
# End Source File
# Begin Source File

SOURCE=.\Shell.cpp
# End Source File
# Begin Source File

SOURCE=.\ShellSM.cpp
# End Source File
# Begin Source File

SOURCE=.\SphereEvaluator.cpp
# End Source File
# Begin Source File

SOURCE=.\SplitSurfaceTool.cpp
# End Source File
# Begin Source File

SOURCE=.\Surface.cpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceOverlapFacet.cpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceOverlapTool.cpp
# End Source File
# Begin Source File

SOURCE=.\SurfaceSM.cpp
# End Source File
# Begin Source File

SOURCE=.\SurfParamTool.cpp
# End Source File
# Begin Source File

SOURCE=.\TDCAGE.cpp
# End Source File
# Begin Source File

SOURCE=.\TDChordal.cpp
# End Source File
# Begin Source File

SOURCE=.\TDParallel.cpp
# End Source File
# Begin Source File

SOURCE=.\TDSplitSurface.cpp
# End Source File
# Begin Source File

SOURCE=.\TDSurfaceOverlap.cpp
# End Source File
# Begin Source File

SOURCE=.\TDUniqueId.cpp
# End Source File
# Begin Source File

SOURCE=.\TopologyBridge.cpp
# End Source File
# Begin Source File

SOURCE=.\TopologyEntity.cpp
# End Source File
# End Group
# Begin Group "facetbool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\facetbool\FBClassify.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBDataUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBImprint.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBIntersect.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBPolyhedron.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBRetriangulate.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\FBTiler.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\IntegerHash.cpp
# End Source File
# Begin Source File

SOURCE=.\facetbool\KdTree.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=.\Readme.txt
# End Source File
# End Target
# End Project
