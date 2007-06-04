/* this ALWAYS GENERATED file contains the definitions for the interfaces */


/* File created by MIDL compiler version 5.01.0164 */
/* at Thu Oct 04 16:54:01 2001
 */
/* Compiler settings for C:\work\gttest\gtfordm.idl:
    Os (OptLev=s), W1, Zp8, env=Win32, ms_ext, c_ext
    error checks: allocation ref bounds_check enum stub_data 
*/
//@@MIDL_FILE_HEADING(  )


/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 440
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif // __RPCNDR_H_VERSION__

#ifndef COM_NO_WINDOWS_H
#include "windows.h"
#include "ole2.h"
#endif /*COM_NO_WINDOWS_H*/

#ifndef __gtfordm_h__
#define __gtfordm_h__

#ifdef __cplusplus
extern "C"{
#endif 

/* Forward Declarations */ 

#ifndef __IDMSurfaceBodies_FWD_DEFINED__
#define __IDMSurfaceBodies_FWD_DEFINED__
typedef interface IDMSurfaceBodies IDMSurfaceBodies;
#endif 	/* __IDMSurfaceBodies_FWD_DEFINED__ */


#ifndef __IDMSurfaceBody_FWD_DEFINED__
#define __IDMSurfaceBody_FWD_DEFINED__
typedef interface IDMSurfaceBody IDMSurfaceBody;
#endif 	/* __IDMSurfaceBody_FWD_DEFINED__ */


#ifndef __IDMShell_FWD_DEFINED__
#define __IDMShell_FWD_DEFINED__
typedef interface IDMShell IDMShell;
#endif 	/* __IDMShell_FWD_DEFINED__ */


#ifndef __IDMFace_FWD_DEFINED__
#define __IDMFace_FWD_DEFINED__
typedef interface IDMFace IDMFace;
#endif 	/* __IDMFace_FWD_DEFINED__ */


#ifndef __IDMLoop_FWD_DEFINED__
#define __IDMLoop_FWD_DEFINED__
typedef interface IDMLoop IDMLoop;
#endif 	/* __IDMLoop_FWD_DEFINED__ */


#ifndef __IDMEdgeUse_FWD_DEFINED__
#define __IDMEdgeUse_FWD_DEFINED__
typedef interface IDMEdgeUse IDMEdgeUse;
#endif 	/* __IDMEdgeUse_FWD_DEFINED__ */


#ifndef __IDMEdge_FWD_DEFINED__
#define __IDMEdge_FWD_DEFINED__
typedef interface IDMEdge IDMEdge;
#endif 	/* __IDMEdge_FWD_DEFINED__ */


#ifndef __IDMVertex_FWD_DEFINED__
#define __IDMVertex_FWD_DEFINED__
typedef interface IDMVertex IDMVertex;
#endif 	/* __IDMVertex_FWD_DEFINED__ */


#ifndef __IDMSurface_FWD_DEFINED__
#define __IDMSurface_FWD_DEFINED__
typedef interface IDMSurface IDMSurface;
#endif 	/* __IDMSurface_FWD_DEFINED__ */


#ifndef __IDMCurve_FWD_DEFINED__
#define __IDMCurve_FWD_DEFINED__
typedef interface IDMCurve IDMCurve;
#endif 	/* __IDMCurve_FWD_DEFINED__ */


#ifndef __IDMCurve2D_FWD_DEFINED__
#define __IDMCurve2D_FWD_DEFINED__
typedef interface IDMCurve2D IDMCurve2D;
#endif 	/* __IDMCurve2D_FWD_DEFINED__ */


#ifndef __IDMCone_FWD_DEFINED__
#define __IDMCone_FWD_DEFINED__
typedef interface IDMCone IDMCone;
#endif 	/* __IDMCone_FWD_DEFINED__ */


#ifndef __IDMCylinder_FWD_DEFINED__
#define __IDMCylinder_FWD_DEFINED__
typedef interface IDMCylinder IDMCylinder;
#endif 	/* __IDMCylinder_FWD_DEFINED__ */


#ifndef __IDMSphere_FWD_DEFINED__
#define __IDMSphere_FWD_DEFINED__
typedef interface IDMSphere IDMSphere;
#endif 	/* __IDMSphere_FWD_DEFINED__ */


#ifndef __IDMTorus_FWD_DEFINED__
#define __IDMTorus_FWD_DEFINED__
typedef interface IDMTorus IDMTorus;
#endif 	/* __IDMTorus_FWD_DEFINED__ */


#ifndef __IDMPlane_FWD_DEFINED__
#define __IDMPlane_FWD_DEFINED__
typedef interface IDMPlane IDMPlane;
#endif 	/* __IDMPlane_FWD_DEFINED__ */


#ifndef __IDMBSplineSurface_FWD_DEFINED__
#define __IDMBSplineSurface_FWD_DEFINED__
typedef interface IDMBSplineSurface IDMBSplineSurface;
#endif 	/* __IDMBSplineSurface_FWD_DEFINED__ */


#ifndef __IDMCircle_FWD_DEFINED__
#define __IDMCircle_FWD_DEFINED__
typedef interface IDMCircle IDMCircle;
#endif 	/* __IDMCircle_FWD_DEFINED__ */


#ifndef __IDMEllipse_FWD_DEFINED__
#define __IDMEllipse_FWD_DEFINED__
typedef interface IDMEllipse IDMEllipse;
#endif 	/* __IDMEllipse_FWD_DEFINED__ */


#ifndef __IDMLine_FWD_DEFINED__
#define __IDMLine_FWD_DEFINED__
typedef interface IDMLine IDMLine;
#endif 	/* __IDMLine_FWD_DEFINED__ */


#ifndef __IDMPolyline_FWD_DEFINED__
#define __IDMPolyline_FWD_DEFINED__
typedef interface IDMPolyline IDMPolyline;
#endif 	/* __IDMPolyline_FWD_DEFINED__ */


#ifndef __IDMBSplineCurve_FWD_DEFINED__
#define __IDMBSplineCurve_FWD_DEFINED__
typedef interface IDMBSplineCurve IDMBSplineCurve;
#endif 	/* __IDMBSplineCurve_FWD_DEFINED__ */


#ifndef __IDMCircle2D_FWD_DEFINED__
#define __IDMCircle2D_FWD_DEFINED__
typedef interface IDMCircle2D IDMCircle2D;
#endif 	/* __IDMCircle2D_FWD_DEFINED__ */


#ifndef __IDMEllipse2D_FWD_DEFINED__
#define __IDMEllipse2D_FWD_DEFINED__
typedef interface IDMEllipse2D IDMEllipse2D;
#endif 	/* __IDMEllipse2D_FWD_DEFINED__ */


#ifndef __IDMLine2D_FWD_DEFINED__
#define __IDMLine2D_FWD_DEFINED__
typedef interface IDMLine2D IDMLine2D;
#endif 	/* __IDMLine2D_FWD_DEFINED__ */


#ifndef __IDMPolyline2D_FWD_DEFINED__
#define __IDMPolyline2D_FWD_DEFINED__
typedef interface IDMPolyline2D IDMPolyline2D;
#endif 	/* __IDMPolyline2D_FWD_DEFINED__ */


#ifndef __IDMBSplineCurve2D_FWD_DEFINED__
#define __IDMBSplineCurve2D_FWD_DEFINED__
typedef interface IDMBSplineCurve2D IDMBSplineCurve2D;
#endif 	/* __IDMBSplineCurve2D_FWD_DEFINED__ */


#ifndef __IEnumDMSurfaceBodies_FWD_DEFINED__
#define __IEnumDMSurfaceBodies_FWD_DEFINED__
typedef interface IEnumDMSurfaceBodies IEnumDMSurfaceBodies;
#endif 	/* __IEnumDMSurfaceBodies_FWD_DEFINED__ */


#ifndef __IEnumDMShells_FWD_DEFINED__
#define __IEnumDMShells_FWD_DEFINED__
typedef interface IEnumDMShells IEnumDMShells;
#endif 	/* __IEnumDMShells_FWD_DEFINED__ */


#ifndef __IEnumDMFaces_FWD_DEFINED__
#define __IEnumDMFaces_FWD_DEFINED__
typedef interface IEnumDMFaces IEnumDMFaces;
#endif 	/* __IEnumDMFaces_FWD_DEFINED__ */


#ifndef __IEnumDMLoops_FWD_DEFINED__
#define __IEnumDMLoops_FWD_DEFINED__
typedef interface IEnumDMLoops IEnumDMLoops;
#endif 	/* __IEnumDMLoops_FWD_DEFINED__ */


#ifndef __IEnumDMEdgeUses_FWD_DEFINED__
#define __IEnumDMEdgeUses_FWD_DEFINED__
typedef interface IEnumDMEdgeUses IEnumDMEdgeUses;
#endif 	/* __IEnumDMEdgeUses_FWD_DEFINED__ */


#ifndef __IEnumDMEdges_FWD_DEFINED__
#define __IEnumDMEdges_FWD_DEFINED__
typedef interface IEnumDMEdges IEnumDMEdges;
#endif 	/* __IEnumDMEdges_FWD_DEFINED__ */


#ifndef __IEnumDMVertices_FWD_DEFINED__
#define __IEnumDMVertices_FWD_DEFINED__
typedef interface IEnumDMVertices IEnumDMVertices;
#endif 	/* __IEnumDMVertices_FWD_DEFINED__ */


#ifndef __IDMAlternateSurfaceBody_FWD_DEFINED__
#define __IDMAlternateSurfaceBody_FWD_DEFINED__
typedef interface IDMAlternateSurfaceBody IDMAlternateSurfaceBody;
#endif 	/* __IDMAlternateSurfaceBody_FWD_DEFINED__ */


#ifndef __IDMReferenceKey_FWD_DEFINED__
#define __IDMReferenceKey_FWD_DEFINED__
typedef interface IDMReferenceKey IDMReferenceKey;
#endif 	/* __IDMReferenceKey_FWD_DEFINED__ */


#ifndef __IDMReference_FWD_DEFINED__
#define __IDMReference_FWD_DEFINED__
typedef interface IDMReference IDMReference;
#endif 	/* __IDMReference_FWD_DEFINED__ */


#ifndef __IEnumDMReferenceKeys_FWD_DEFINED__
#define __IEnumDMReferenceKeys_FWD_DEFINED__
typedef interface IEnumDMReferenceKeys IEnumDMReferenceKeys;
#endif 	/* __IEnumDMReferenceKeys_FWD_DEFINED__ */


#ifndef __IDMGeometricLocate_FWD_DEFINED__
#define __IDMGeometricLocate_FWD_DEFINED__
typedef interface IDMGeometricLocate IDMGeometricLocate;
#endif 	/* __IDMGeometricLocate_FWD_DEFINED__ */


/* header files for imported files */
#include "unknwn.h"
#include "objidl.h"

void __RPC_FAR * __RPC_USER MIDL_user_allocate(size_t);
void __RPC_USER MIDL_user_free( void __RPC_FAR * ); 

/* interface __MIDL_itf_gtfordm_0000 */
/* [local] */ 

#ifndef __WINDOWS_H__
#include <windows.h>
#endif
#ifndef __OBJBASE_H__
#include <objbase.h>
#endif
DEFINE_GUID(IID_IDMSurfaceBodies, 0x0002D208, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMSurfaceBody, 0x0002D209, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMShell, 0x0002D20A, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMFace, 0x0002D20B, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMLoop, 0x0002D20C, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMEdgeUse, 0x0002D20D, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMEdge, 0x0002D20E, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMVertex, 0x0002D20F, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMSurface, 0x0002D210, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCurve, 0x0002D211, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCurve2D, 0x0002D212, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCone, 0x0002D213, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCylinder, 0x0002D214, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMSphere, 0x0002D215, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMTorus, 0x0002D216, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMPlane, 0x0002D217, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMBSplineSurface, 0x0002D218, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCircle, 0x0002D219, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMEllipse, 0x0002D21A, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMLine, 0x0002D21B, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMPolyLine, 0x0002D21C, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMBSplineCurve, 0x0002D21D, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMCircle2D, 0x0002D21E, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMEllipse2D, 0x0002D21F, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMLine2D, 0x0002D220, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMPolyLine2D, 0x0002D221, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMBSplineCurve2D, 0x0002D222, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMSurfaceBodies, 0x0002D223, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMShells, 0x0002D224, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMFaces, 0x0002D225, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMLoops, 0x0002D226, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMEdgeUses, 0x0002D227, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMEdges, 0x0002D228, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMVertices, 0x0002D229, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMAlternateSurfaceBody, 0x0002D22A, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
typedef /* [unique] */ IDMSurfaceBodies __RPC_FAR *LPDMSURFACEBODIES;

typedef /* [unique] */ IDMSurfaceBody __RPC_FAR *LPDMSURFACEBODY;

typedef /* [unique] */ IDMShell __RPC_FAR *LPDMSHELL;

typedef /* [unique] */ IDMFace __RPC_FAR *LPDMFACE;

typedef /* [unique] */ IDMLoop __RPC_FAR *LPDMLOOP;

typedef /* [unique] */ IDMEdgeUse __RPC_FAR *LPDMEDGEUSE;

typedef /* [unique] */ IDMEdge __RPC_FAR *LPDMEDGE;

typedef /* [unique] */ IDMVertex __RPC_FAR *LPDMVERTEX;

typedef /* [unique] */ IDMSurface __RPC_FAR *LPDMSURFACE;

typedef /* [unique] */ IDMCurve __RPC_FAR *LPDMCURVE;

typedef /* [unique] */ IDMCurve2D __RPC_FAR *LPDMCURVE2D;

typedef /* [unique] */ IDMCone __RPC_FAR *LPDMCONE;

typedef /* [unique] */ IDMCylinder __RPC_FAR *LPDMCYLINDER;

typedef /* [unique] */ IDMSphere __RPC_FAR *LPDMSPHERE;

typedef /* [unique] */ IDMTorus __RPC_FAR *LPDMTORUS;

typedef /* [unique] */ IDMPlane __RPC_FAR *LPDMPLANE;

typedef /* [unique] */ IDMBSplineSurface __RPC_FAR *LPDMBSPLINESURFACE;

typedef /* [unique] */ IDMCircle __RPC_FAR *LPDMCIRCLE;

typedef /* [unique] */ IDMEllipse __RPC_FAR *LPDMELLIPSE;

typedef /* [unique] */ IDMLine __RPC_FAR *LPDMLINE;

typedef /* [unique] */ IDMPolyline __RPC_FAR *LPDMPOLYLINE;

typedef /* [unique] */ IDMBSplineCurve __RPC_FAR *LPDMBSPLINECURVE;

typedef /* [unique] */ IDMCircle2D __RPC_FAR *LPDMCIRCLE2D;

typedef /* [unique] */ IDMEllipse2D __RPC_FAR *LPDMELLIPSE2D;

typedef /* [unique] */ IDMLine2D __RPC_FAR *LPDMLINE2D;

typedef /* [unique] */ IDMPolyline2D __RPC_FAR *LPDMPOLYLINE2D;

typedef /* [unique] */ IDMBSplineCurve2D __RPC_FAR *LPDMBSPLINECURVE2D;

typedef /* [unique] */ IEnumDMSurfaceBodies __RPC_FAR *LPENUM_DMSURFACEBODIES;

typedef /* [unique] */ IEnumDMShells __RPC_FAR *LPENUM_DMSHELLS;

typedef /* [unique] */ IEnumDMFaces __RPC_FAR *LPENUM_DMFACES;

typedef /* [unique] */ IEnumDMLoops __RPC_FAR *LPENUM_DMLOOPS;

typedef /* [unique] */ IEnumDMEdgeUses __RPC_FAR *LPENUM_DMEDGEUSES;

typedef /* [unique] */ IEnumDMEdges __RPC_FAR *LPENUM_DMEDGES;

typedef /* [unique] */ IEnumDMVertices __RPC_FAR *LPENUM_DMVERTICES;

typedef /* [unique] */ IDMAlternateSurfaceBody __RPC_FAR *LPDMALTERNATESURFACEBODY;

typedef 
enum tagDMSFGEOMETRYFORM
    {	SFGEOMETRYFORM_CLOSEDUVLOOPS	= 0x1,
	SFGEOMETRYFORM_NOT_CLOSEDUVLOOPS	= 0x2,
	SFGEOMETRYFORM_NURBS	= 0x4,
	SFGEOMETRYFORM_NOT_NURBS	= 0x8
    }	DMSFGEOMETRYFORM;

typedef 
enum tagDMCVGEOMETRYFORM
    {	CVGEOMETRYFORM_NURBS	= 0x1,
	CVGEOMETRYFORM_NOT_NURBS	= 0x2
    }	DMCVGEOMETRYFORM;



extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0000_v0_0_c_ifspec;
extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0000_v0_0_s_ifspec;

#ifndef __IDMSurfaceBodies_INTERFACE_DEFINED__
#define __IDMSurfaceBodies_INTERFACE_DEFINED__

/* interface IDMSurfaceBodies */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMSurfaceBodies;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D208-0000-0000-C000-000000000046")
    IDMSurfaceBodies : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumSurfaceBodies( 
            /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMSurfaceBodiesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMSurfaceBodies __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMSurfaceBodies __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMSurfaceBodies __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumSurfaceBodies )( 
            IDMSurfaceBodies __RPC_FAR * This,
            /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies);
        
        END_INTERFACE
    } IDMSurfaceBodiesVtbl;

    interface IDMSurfaceBodies
    {
        CONST_VTBL struct IDMSurfaceBodiesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMSurfaceBodies_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMSurfaceBodies_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMSurfaceBodies_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMSurfaceBodies_EnumSurfaceBodies(This,pEnumSurfaceBodies)	\
    (This)->lpVtbl -> EnumSurfaceBodies(This,pEnumSurfaceBodies)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMSurfaceBodies_EnumSurfaceBodies_Proxy( 
    IDMSurfaceBodies __RPC_FAR * This,
    /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies);


void __RPC_STUB IDMSurfaceBodies_EnumSurfaceBodies_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMSurfaceBodies_INTERFACE_DEFINED__ */


#ifndef __IDMSurfaceBody_INTERFACE_DEFINED__
#define __IDMSurfaceBody_INTERFACE_DEFINED__

/* interface IDMSurfaceBody */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMSurfaceBody;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D209-0000-0000-C000-000000000046")
    IDMSurfaceBody : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumShells( 
            /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumFaces( 
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumEdges( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetDocument( 
            /* [out] */ LPDMSURFACEBODIES __RPC_FAR *pDoc) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsSolid( 
            /* [out] */ boolean __RPC_FAR *pIsSolid) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetVolume( 
            /* [out] */ double __RPC_FAR *pVolume) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetGeometryForm( 
            /* [out] */ DWORD __RPC_FAR *pForm) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMSurfaceBodyVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMSurfaceBody __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMSurfaceBody __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumShells )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumFaces )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdges )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetDocument )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ LPDMSURFACEBODIES __RPC_FAR *pDoc);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsSolid )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsSolid);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetVolume )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *pVolume);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetGeometryForm )( 
            IDMSurfaceBody __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *pForm);
        
        END_INTERFACE
    } IDMSurfaceBodyVtbl;

    interface IDMSurfaceBody
    {
        CONST_VTBL struct IDMSurfaceBodyVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMSurfaceBody_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMSurfaceBody_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMSurfaceBody_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMSurfaceBody_EnumShells(This,pEnumShells)	\
    (This)->lpVtbl -> EnumShells(This,pEnumShells)

#define IDMSurfaceBody_EnumFaces(This,pEnumFaces)	\
    (This)->lpVtbl -> EnumFaces(This,pEnumFaces)

#define IDMSurfaceBody_EnumEdges(This,pEnumEdges)	\
    (This)->lpVtbl -> EnumEdges(This,pEnumEdges)

#define IDMSurfaceBody_GetDocument(This,pDoc)	\
    (This)->lpVtbl -> GetDocument(This,pDoc)

#define IDMSurfaceBody_IsSolid(This,pIsSolid)	\
    (This)->lpVtbl -> IsSolid(This,pIsSolid)

#define IDMSurfaceBody_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#define IDMSurfaceBody_GetVolume(This,pVolume)	\
    (This)->lpVtbl -> GetVolume(This,pVolume)

#define IDMSurfaceBody_GetGeometryForm(This,pForm)	\
    (This)->lpVtbl -> GetGeometryForm(This,pForm)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMSurfaceBody_EnumShells_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells);


void __RPC_STUB IDMSurfaceBody_EnumShells_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_EnumFaces_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);


void __RPC_STUB IDMSurfaceBody_EnumFaces_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_EnumEdges_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IDMSurfaceBody_EnumEdges_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_GetDocument_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ LPDMSURFACEBODIES __RPC_FAR *pDoc);


void __RPC_STUB IDMSurfaceBody_GetDocument_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_IsSolid_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsSolid);


void __RPC_STUB IDMSurfaceBody_IsSolid_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_GetRangeBox_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);


void __RPC_STUB IDMSurfaceBody_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_GetVolume_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *pVolume);


void __RPC_STUB IDMSurfaceBody_GetVolume_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurfaceBody_GetGeometryForm_Proxy( 
    IDMSurfaceBody __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *pForm);


void __RPC_STUB IDMSurfaceBody_GetGeometryForm_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMSurfaceBody_INTERFACE_DEFINED__ */


#ifndef __IDMShell_INTERFACE_DEFINED__
#define __IDMShell_INTERFACE_DEFINED__

/* interface IDMShell */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMShell;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20A-0000-0000-C000-000000000046")
    IDMShell : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumFaces( 
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumEdges( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetBody( 
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsClosed( 
            /* [out] */ boolean __RPC_FAR *pIsClosed) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsVoid( 
            /* [out] */ boolean __RPC_FAR *pIsVoid) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsPointInside( 
            /* [in] */ double __RPC_FAR point[ 3 ],
            /* [out] */ boolean __RPC_FAR *pIsPtInside) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetVolume( 
            /* [out] */ double __RPC_FAR *pVolume) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMShellVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMShell __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMShell __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMShell __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumFaces )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdges )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBody )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsClosed )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsClosed);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsVoid )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsVoid);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsPointInside )( 
            IDMShell __RPC_FAR * This,
            /* [in] */ double __RPC_FAR point[ 3 ],
            /* [out] */ boolean __RPC_FAR *pIsPtInside);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetVolume )( 
            IDMShell __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *pVolume);
        
        END_INTERFACE
    } IDMShellVtbl;

    interface IDMShell
    {
        CONST_VTBL struct IDMShellVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMShell_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMShell_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMShell_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMShell_EnumFaces(This,pEnumFaces)	\
    (This)->lpVtbl -> EnumFaces(This,pEnumFaces)

#define IDMShell_EnumEdges(This,pEnumEdges)	\
    (This)->lpVtbl -> EnumEdges(This,pEnumEdges)

#define IDMShell_GetBody(This,pBody)	\
    (This)->lpVtbl -> GetBody(This,pBody)

#define IDMShell_IsClosed(This,pIsClosed)	\
    (This)->lpVtbl -> IsClosed(This,pIsClosed)

#define IDMShell_IsVoid(This,pIsVoid)	\
    (This)->lpVtbl -> IsVoid(This,pIsVoid)

#define IDMShell_IsPointInside(This,point,pIsPtInside)	\
    (This)->lpVtbl -> IsPointInside(This,point,pIsPtInside)

#define IDMShell_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#define IDMShell_GetVolume(This,pVolume)	\
    (This)->lpVtbl -> GetVolume(This,pVolume)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMShell_EnumFaces_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);


void __RPC_STUB IDMShell_EnumFaces_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_EnumEdges_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IDMShell_EnumEdges_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_GetBody_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody);


void __RPC_STUB IDMShell_GetBody_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_IsClosed_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsClosed);


void __RPC_STUB IDMShell_IsClosed_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_IsVoid_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsVoid);


void __RPC_STUB IDMShell_IsVoid_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_IsPointInside_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [in] */ double __RPC_FAR point[ 3 ],
    /* [out] */ boolean __RPC_FAR *pIsPtInside);


void __RPC_STUB IDMShell_IsPointInside_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_GetRangeBox_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);


void __RPC_STUB IDMShell_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMShell_GetVolume_Proxy( 
    IDMShell __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *pVolume);


void __RPC_STUB IDMShell_GetVolume_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMShell_INTERFACE_DEFINED__ */


#ifndef __IDMFace_INTERFACE_DEFINED__
#define __IDMFace_INTERFACE_DEFINED__

/* interface IDMFace */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMFace;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20B-0000-0000-C000-000000000046")
    IDMFace : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumLoops( 
            /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumEdges( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumVertices( 
            /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetShell( 
            /* [out] */ LPDMSHELL __RPC_FAR *pShell) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetBody( 
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsParamReversed( 
            /* [out] */ boolean __RPC_FAR *pIsReversed) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMFaceVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMFace __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMFace __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMFace __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumLoops )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdges )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumVertices )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetShell )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ LPDMSHELL __RPC_FAR *pShell);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBody )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsParamReversed )( 
            IDMFace __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsReversed);
        
        END_INTERFACE
    } IDMFaceVtbl;

    interface IDMFace
    {
        CONST_VTBL struct IDMFaceVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMFace_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMFace_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMFace_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMFace_EnumLoops(This,pEnumLoops)	\
    (This)->lpVtbl -> EnumLoops(This,pEnumLoops)

#define IDMFace_EnumEdges(This,pEnumEdges)	\
    (This)->lpVtbl -> EnumEdges(This,pEnumEdges)

#define IDMFace_EnumVertices(This,pEnumVertices)	\
    (This)->lpVtbl -> EnumVertices(This,pEnumVertices)

#define IDMFace_GetShell(This,pShell)	\
    (This)->lpVtbl -> GetShell(This,pShell)

#define IDMFace_GetBody(This,pBody)	\
    (This)->lpVtbl -> GetBody(This,pBody)

#define IDMFace_IsParamReversed(This,pIsReversed)	\
    (This)->lpVtbl -> IsParamReversed(This,pIsReversed)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMFace_EnumLoops_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops);


void __RPC_STUB IDMFace_EnumLoops_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMFace_EnumEdges_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IDMFace_EnumEdges_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMFace_EnumVertices_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices);


void __RPC_STUB IDMFace_EnumVertices_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMFace_GetShell_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ LPDMSHELL __RPC_FAR *pShell);


void __RPC_STUB IDMFace_GetShell_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMFace_GetBody_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ LPDMSURFACEBODY __RPC_FAR *pBody);


void __RPC_STUB IDMFace_GetBody_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMFace_IsParamReversed_Proxy( 
    IDMFace __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsReversed);


void __RPC_STUB IDMFace_IsParamReversed_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMFace_INTERFACE_DEFINED__ */


#ifndef __IDMLoop_INTERFACE_DEFINED__
#define __IDMLoop_INTERFACE_DEFINED__

/* interface IDMLoop */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMLoop;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20C-0000-0000-C000-000000000046")
    IDMLoop : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumEdgeUses( 
            /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumEdges( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetFace( 
            /* [out] */ LPDMFACE __RPC_FAR *pFace) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsOuterLoop( 
            /* [out] */ boolean __RPC_FAR *pIsOuter) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMLoopVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMLoop __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMLoop __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMLoop __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdgeUses )( 
            IDMLoop __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdges )( 
            IDMLoop __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetFace )( 
            IDMLoop __RPC_FAR * This,
            /* [out] */ LPDMFACE __RPC_FAR *pFace);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsOuterLoop )( 
            IDMLoop __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsOuter);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMLoop __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);
        
        END_INTERFACE
    } IDMLoopVtbl;

    interface IDMLoop
    {
        CONST_VTBL struct IDMLoopVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMLoop_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMLoop_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMLoop_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMLoop_EnumEdgeUses(This,pEnumEdgeUses)	\
    (This)->lpVtbl -> EnumEdgeUses(This,pEnumEdgeUses)

#define IDMLoop_EnumEdges(This,pEnumEdges)	\
    (This)->lpVtbl -> EnumEdges(This,pEnumEdges)

#define IDMLoop_GetFace(This,pFace)	\
    (This)->lpVtbl -> GetFace(This,pFace)

#define IDMLoop_IsOuterLoop(This,pIsOuter)	\
    (This)->lpVtbl -> IsOuterLoop(This,pIsOuter)

#define IDMLoop_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMLoop_EnumEdgeUses_Proxy( 
    IDMLoop __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses);


void __RPC_STUB IDMLoop_EnumEdgeUses_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMLoop_EnumEdges_Proxy( 
    IDMLoop __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IDMLoop_EnumEdges_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMLoop_GetFace_Proxy( 
    IDMLoop __RPC_FAR * This,
    /* [out] */ LPDMFACE __RPC_FAR *pFace);


void __RPC_STUB IDMLoop_GetFace_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMLoop_IsOuterLoop_Proxy( 
    IDMLoop __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsOuter);


void __RPC_STUB IDMLoop_IsOuterLoop_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMLoop_GetRangeBox_Proxy( 
    IDMLoop __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);


void __RPC_STUB IDMLoop_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMLoop_INTERFACE_DEFINED__ */


#ifndef __IDMEdgeUse_INTERFACE_DEFINED__
#define __IDMEdgeUse_INTERFACE_DEFINED__

/* interface IDMEdgeUse */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMEdgeUse;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20D-0000-0000-C000-000000000046")
    IDMEdgeUse : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetEdge( 
            /* [out] */ LPDMEDGE __RPC_FAR *pEdge) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetLoop( 
            /* [out] */ LPDMLOOP __RPC_FAR *pLoop) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPartner( 
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetNext( 
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPrevious( 
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsOpposedToEdge( 
            /* [out] */ boolean __RPC_FAR *pIsReversed) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsParamReversed( 
            /* [out] */ boolean __RPC_FAR *pIsParamReversed) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMEdgeUseVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMEdgeUse __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMEdgeUse __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEdge )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ LPDMEDGE __RPC_FAR *pEdge);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetLoop )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ LPDMLOOP __RPC_FAR *pLoop);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPartner )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetNext )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPrevious )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsOpposedToEdge )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsReversed);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsParamReversed )( 
            IDMEdgeUse __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsParamReversed);
        
        END_INTERFACE
    } IDMEdgeUseVtbl;

    interface IDMEdgeUse
    {
        CONST_VTBL struct IDMEdgeUseVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMEdgeUse_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMEdgeUse_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMEdgeUse_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMEdgeUse_GetEdge(This,pEdge)	\
    (This)->lpVtbl -> GetEdge(This,pEdge)

#define IDMEdgeUse_GetLoop(This,pLoop)	\
    (This)->lpVtbl -> GetLoop(This,pLoop)

#define IDMEdgeUse_GetPartner(This,pEdgeUse)	\
    (This)->lpVtbl -> GetPartner(This,pEdgeUse)

#define IDMEdgeUse_GetNext(This,pEdgeUse)	\
    (This)->lpVtbl -> GetNext(This,pEdgeUse)

#define IDMEdgeUse_GetPrevious(This,pEdgeUse)	\
    (This)->lpVtbl -> GetPrevious(This,pEdgeUse)

#define IDMEdgeUse_IsOpposedToEdge(This,pIsReversed)	\
    (This)->lpVtbl -> IsOpposedToEdge(This,pIsReversed)

#define IDMEdgeUse_IsParamReversed(This,pIsParamReversed)	\
    (This)->lpVtbl -> IsParamReversed(This,pIsParamReversed)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMEdgeUse_GetEdge_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ LPDMEDGE __RPC_FAR *pEdge);


void __RPC_STUB IDMEdgeUse_GetEdge_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_GetLoop_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ LPDMLOOP __RPC_FAR *pLoop);


void __RPC_STUB IDMEdgeUse_GetLoop_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_GetPartner_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);


void __RPC_STUB IDMEdgeUse_GetPartner_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_GetNext_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);


void __RPC_STUB IDMEdgeUse_GetNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_GetPrevious_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse);


void __RPC_STUB IDMEdgeUse_GetPrevious_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_IsOpposedToEdge_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsReversed);


void __RPC_STUB IDMEdgeUse_IsOpposedToEdge_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdgeUse_IsParamReversed_Proxy( 
    IDMEdgeUse __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsParamReversed);


void __RPC_STUB IDMEdgeUse_IsParamReversed_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMEdgeUse_INTERFACE_DEFINED__ */


#ifndef __IDMEdge_INTERFACE_DEFINED__
#define __IDMEdge_INTERFACE_DEFINED__

/* interface IDMEdge */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMEdge;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20E-0000-0000-C000-000000000046")
    IDMEdge : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetStartVertex( 
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetEndVertex( 
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetEdgeUses( 
            /* [out] */ ULONG __RPC_FAR *pNumEdgeUses,
            /* [out] */ LPDMEDGEUSE __RPC_FAR pEdgeUses[ 2 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetFaces( 
            /* [out] */ ULONG __RPC_FAR *pNumFaces,
            /* [out] */ LPDMFACE __RPC_FAR pFaces[ 2 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsParamReversed( 
            /* [out] */ boolean __RPC_FAR *pIsParamReversed) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMEdgeVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMEdge __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMEdge __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMEdge __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetStartVertex )( 
            IDMEdge __RPC_FAR * This,
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEndVertex )( 
            IDMEdge __RPC_FAR * This,
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEdgeUses )( 
            IDMEdge __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pNumEdgeUses,
            /* [out] */ LPDMEDGEUSE __RPC_FAR pEdgeUses[ 2 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetFaces )( 
            IDMEdge __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pNumFaces,
            /* [out] */ LPDMFACE __RPC_FAR pFaces[ 2 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsParamReversed )( 
            IDMEdge __RPC_FAR * This,
            /* [out] */ boolean __RPC_FAR *pIsParamReversed);
        
        END_INTERFACE
    } IDMEdgeVtbl;

    interface IDMEdge
    {
        CONST_VTBL struct IDMEdgeVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMEdge_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMEdge_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMEdge_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMEdge_GetStartVertex(This,pVertex)	\
    (This)->lpVtbl -> GetStartVertex(This,pVertex)

#define IDMEdge_GetEndVertex(This,pVertex)	\
    (This)->lpVtbl -> GetEndVertex(This,pVertex)

#define IDMEdge_GetEdgeUses(This,pNumEdgeUses,pEdgeUses)	\
    (This)->lpVtbl -> GetEdgeUses(This,pNumEdgeUses,pEdgeUses)

#define IDMEdge_GetFaces(This,pNumFaces,pFaces)	\
    (This)->lpVtbl -> GetFaces(This,pNumFaces,pFaces)

#define IDMEdge_IsParamReversed(This,pIsParamReversed)	\
    (This)->lpVtbl -> IsParamReversed(This,pIsParamReversed)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMEdge_GetStartVertex_Proxy( 
    IDMEdge __RPC_FAR * This,
    /* [out] */ LPDMVERTEX __RPC_FAR *pVertex);


void __RPC_STUB IDMEdge_GetStartVertex_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdge_GetEndVertex_Proxy( 
    IDMEdge __RPC_FAR * This,
    /* [out] */ LPDMVERTEX __RPC_FAR *pVertex);


void __RPC_STUB IDMEdge_GetEndVertex_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdge_GetEdgeUses_Proxy( 
    IDMEdge __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pNumEdgeUses,
    /* [out] */ LPDMEDGEUSE __RPC_FAR pEdgeUses[ 2 ]);


void __RPC_STUB IDMEdge_GetEdgeUses_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdge_GetFaces_Proxy( 
    IDMEdge __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pNumFaces,
    /* [out] */ LPDMFACE __RPC_FAR pFaces[ 2 ]);


void __RPC_STUB IDMEdge_GetFaces_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMEdge_IsParamReversed_Proxy( 
    IDMEdge __RPC_FAR * This,
    /* [out] */ boolean __RPC_FAR *pIsParamReversed);


void __RPC_STUB IDMEdge_IsParamReversed_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMEdge_INTERFACE_DEFINED__ */


#ifndef __IDMVertex_INTERFACE_DEFINED__
#define __IDMVertex_INTERFACE_DEFINED__

/* interface IDMVertex */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMVertex;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D20F-0000-0000-C000-000000000046")
    IDMVertex : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE EnumEdges( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE EnumFaces( 
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPoint( 
            /* [out] */ double __RPC_FAR point[ 3 ]) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMVertexVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMVertex __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMVertex __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMVertex __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumEdges )( 
            IDMVertex __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *EnumFaces )( 
            IDMVertex __RPC_FAR * This,
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPoint )( 
            IDMVertex __RPC_FAR * This,
            /* [out] */ double __RPC_FAR point[ 3 ]);
        
        END_INTERFACE
    } IDMVertexVtbl;

    interface IDMVertex
    {
        CONST_VTBL struct IDMVertexVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMVertex_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMVertex_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMVertex_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMVertex_EnumEdges(This,pEnumEdges)	\
    (This)->lpVtbl -> EnumEdges(This,pEnumEdges)

#define IDMVertex_EnumFaces(This,pEnumFaces)	\
    (This)->lpVtbl -> EnumFaces(This,pEnumFaces)

#define IDMVertex_GetPoint(This,point)	\
    (This)->lpVtbl -> GetPoint(This,point)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMVertex_EnumEdges_Proxy( 
    IDMVertex __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IDMVertex_EnumEdges_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMVertex_EnumFaces_Proxy( 
    IDMVertex __RPC_FAR * This,
    /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);


void __RPC_STUB IDMVertex_EnumFaces_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMVertex_GetPoint_Proxy( 
    IDMVertex __RPC_FAR * This,
    /* [out] */ double __RPC_FAR point[ 3 ]);


void __RPC_STUB IDMVertex_GetPoint_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMVertex_INTERFACE_DEFINED__ */


#ifndef __IDMSurface_INTERFACE_DEFINED__
#define __IDMSurface_INTERFACE_DEFINED__

/* interface IDMSurface */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMSurface;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D210-0000-0000-C000-000000000046")
    IDMSurface : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetSurfaceType( 
            /* [out] */ IID __RPC_FAR *pIID) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamRangeRect( 
            /* [out] */ double __RPC_FAR pMinParam[ 2 ],
            /* [out] */ double __RPC_FAR pMaxParam[ 2 ]) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetParamAtPoint( 
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPointAtParam( 
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetNormal( 
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pNormals) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetTangents( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pUTangents,
            /* [out] */ double __RPC_FAR *pVTangents) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetCurvatures( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pMaxTangents,
            /* [out] */ double __RPC_FAR *pMaxCurvatures,
            /* [out] */ double __RPC_FAR *pMinCurvatures) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetDerivatives( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pUPartials,
            /* [out] */ double __RPC_FAR *pVPartials,
            /* [out] */ double __RPC_FAR *pUUPartials,
            /* [out] */ double __RPC_FAR *pUVPartials,
            /* [out] */ double __RPC_FAR *pVVPartials,
            /* [out] */ double __RPC_FAR *pUUUPartials,
            /* [out] */ double __RPC_FAR *pVVVPartials) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE IsParamOnFace( 
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ boolean __RPC_FAR *pIsOnFace) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetArea( 
            /* [out] */ double __RPC_FAR *pArea) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetContinuity( 
            /* [out] */ DWORD __RPC_FAR *nLevel) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamAnomaly( 
            /* [out] */ double __RPC_FAR pPeriodicityU[ 2 ],
            /* [out] */ double __RPC_FAR pPeriodicityV[ 2 ],
            /* [out] */ ULONG __RPC_FAR *pnEndSingularityU,
            /* [out] */ double __RPC_FAR pSingularityU[ 2 ],
            /* [out] */ ULONG __RPC_FAR *pnEndSingularityV,
            /* [out] */ double __RPC_FAR pSingularityV[ 2 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetGeometryForm( 
            /* [out] */ DWORD __RPC_FAR *pForm) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMSurfaceVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMSurface __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMSurface __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetSurfaceType )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ IID __RPC_FAR *pIID);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamRangeRect )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinParam[ 2 ],
            /* [out] */ double __RPC_FAR pMaxParam[ 2 ]);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAtPoint )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPointAtParam )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetNormal )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pNormals);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTangents )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pUTangents,
            /* [out] */ double __RPC_FAR *pVTangents);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCurvatures )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pMaxTangents,
            /* [out] */ double __RPC_FAR *pMaxCurvatures,
            /* [out] */ double __RPC_FAR *pMinCurvatures);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetDerivatives )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pUPartials,
            /* [out] */ double __RPC_FAR *pVPartials,
            /* [out] */ double __RPC_FAR *pUUPartials,
            /* [out] */ double __RPC_FAR *pUVPartials,
            /* [out] */ double __RPC_FAR *pVVPartials,
            /* [out] */ double __RPC_FAR *pUUUPartials,
            /* [out] */ double __RPC_FAR *pVVVPartials);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *IsParamOnFace )( 
            IDMSurface __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ boolean __RPC_FAR *pIsOnFace);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetArea )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *pArea);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetContinuity )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *nLevel);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAnomaly )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pPeriodicityU[ 2 ],
            /* [out] */ double __RPC_FAR pPeriodicityV[ 2 ],
            /* [out] */ ULONG __RPC_FAR *pnEndSingularityU,
            /* [out] */ double __RPC_FAR pSingularityU[ 2 ],
            /* [out] */ ULONG __RPC_FAR *pnEndSingularityV,
            /* [out] */ double __RPC_FAR pSingularityV[ 2 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetGeometryForm )( 
            IDMSurface __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *pForm);
        
        END_INTERFACE
    } IDMSurfaceVtbl;

    interface IDMSurface
    {
        CONST_VTBL struct IDMSurfaceVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMSurface_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMSurface_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMSurface_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMSurface_GetSurfaceType(This,pIID)	\
    (This)->lpVtbl -> GetSurfaceType(This,pIID)

#define IDMSurface_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#define IDMSurface_GetParamRangeRect(This,pMinParam,pMaxParam)	\
    (This)->lpVtbl -> GetParamRangeRect(This,pMinParam,pMaxParam)

#define IDMSurface_GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)	\
    (This)->lpVtbl -> GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)

#define IDMSurface_GetPointAtParam(This,nParams,pParams,pPoints)	\
    (This)->lpVtbl -> GetPointAtParam(This,nParams,pParams,pPoints)

#define IDMSurface_GetNormal(This,nParams,pParams,pNormals)	\
    (This)->lpVtbl -> GetNormal(This,nParams,pParams,pNormals)

#define IDMSurface_GetTangents(This,nParams,pParams,pUTangents,pVTangents)	\
    (This)->lpVtbl -> GetTangents(This,nParams,pParams,pUTangents,pVTangents)

#define IDMSurface_GetCurvatures(This,nParams,pParams,pMaxTangents,pMaxCurvatures,pMinCurvatures)	\
    (This)->lpVtbl -> GetCurvatures(This,nParams,pParams,pMaxTangents,pMaxCurvatures,pMinCurvatures)

#define IDMSurface_GetDerivatives(This,nParams,pParams,pUPartials,pVPartials,pUUPartials,pUVPartials,pVVPartials,pUUUPartials,pVVVPartials)	\
    (This)->lpVtbl -> GetDerivatives(This,nParams,pParams,pUPartials,pVPartials,pUUPartials,pUVPartials,pVVPartials,pUUUPartials,pVVVPartials)

#define IDMSurface_IsParamOnFace(This,nParams,pParams,pIsOnFace)	\
    (This)->lpVtbl -> IsParamOnFace(This,nParams,pParams,pIsOnFace)

#define IDMSurface_GetArea(This,pArea)	\
    (This)->lpVtbl -> GetArea(This,pArea)

#define IDMSurface_GetContinuity(This,nLevel)	\
    (This)->lpVtbl -> GetContinuity(This,nLevel)

#define IDMSurface_GetParamAnomaly(This,pPeriodicityU,pPeriodicityV,pnEndSingularityU,pSingularityU,pnEndSingularityV,pSingularityV)	\
    (This)->lpVtbl -> GetParamAnomaly(This,pPeriodicityU,pPeriodicityV,pnEndSingularityU,pSingularityU,pnEndSingularityV,pSingularityV)

#define IDMSurface_GetGeometryForm(This,pForm)	\
    (This)->lpVtbl -> GetGeometryForm(This,pForm)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMSurface_GetSurfaceType_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ IID __RPC_FAR *pIID);


void __RPC_STUB IDMSurface_GetSurfaceType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetRangeBox_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);


void __RPC_STUB IDMSurface_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetParamRangeRect_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinParam[ 2 ],
    /* [out] */ double __RPC_FAR pMaxParam[ 2 ]);


void __RPC_STUB IDMSurface_GetParamRangeRect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_RemoteGetParamAtPoint_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);


void __RPC_STUB IDMSurface_RemoteGetParamAtPoint_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetPointAtParam_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in][size_is] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pPoints);


void __RPC_STUB IDMSurface_GetPointAtParam_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetNormal_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in][size_is] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pNormals);


void __RPC_STUB IDMSurface_GetNormal_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_RemoteGetTangents_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nUTangents,
    /* [size_is][out] */ double __RPC_FAR *pUTangents,
    /* [in] */ ULONG nVTangents,
    /* [size_is][out] */ double __RPC_FAR *pVTangents);


void __RPC_STUB IDMSurface_RemoteGetTangents_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_RemoteGetCurvatures_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nMaxTangents,
    /* [size_is][out] */ double __RPC_FAR *pMaxTangents,
    /* [in] */ ULONG nMaxCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pMaxCurvatures,
    /* [in] */ ULONG nMinCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pMinCurvatures);


void __RPC_STUB IDMSurface_RemoteGetCurvatures_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_RemoteGetDerivatives_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUPartials,
    /* [in] */ ULONG nVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVPartials,
    /* [in] */ ULONG nUUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUUPartials,
    /* [in] */ ULONG nUVPartials,
    /* [size_is][out] */ double __RPC_FAR *pUVPartials,
    /* [in] */ ULONG nVVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVVPartials,
    /* [in] */ ULONG nUUUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUUUPartials,
    /* [in] */ ULONG nVVVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVVVPartials);


void __RPC_STUB IDMSurface_RemoteGetDerivatives_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_IsParamOnFace_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in][size_is] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ boolean __RPC_FAR *pIsOnFace);


void __RPC_STUB IDMSurface_IsParamOnFace_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetArea_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *pArea);


void __RPC_STUB IDMSurface_GetArea_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetContinuity_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *nLevel);


void __RPC_STUB IDMSurface_GetContinuity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetParamAnomaly_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pPeriodicityU[ 2 ],
    /* [out] */ double __RPC_FAR pPeriodicityV[ 2 ],
    /* [out] */ ULONG __RPC_FAR *pnEndSingularityU,
    /* [out] */ double __RPC_FAR pSingularityU[ 2 ],
    /* [out] */ ULONG __RPC_FAR *pnEndSingularityV,
    /* [out] */ double __RPC_FAR pSingularityV[ 2 ]);


void __RPC_STUB IDMSurface_GetParamAnomaly_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMSurface_GetGeometryForm_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *pForm);


void __RPC_STUB IDMSurface_GetGeometryForm_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMSurface_INTERFACE_DEFINED__ */


#ifndef __IDMCurve_INTERFACE_DEFINED__
#define __IDMCurve_INTERFACE_DEFINED__

/* interface IDMCurve */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCurve;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D211-0000-0000-C000-000000000046")
    IDMCurve : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetCurveType( 
            /* [out] */ IID __RPC_FAR *pIID) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetEndPoints( 
            /* [out] */ double __RPC_FAR pStartPoint[ 3 ],
            /* [out] */ double __RPC_FAR pEndPoint[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamExtents( 
            /* [out] */ double __RPC_FAR *pMinParam,
            /* [out] */ double __RPC_FAR *pMaxParam) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetParamAtPoint( 
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPointAtParam( 
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetTangent( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pTangents) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetCurvature( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pDirections,
            /* [out] */ double __RPC_FAR *pCurvatures) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetDerivatives( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pFirstDerivs,
            /* [out] */ double __RPC_FAR *pSecondDerivs,
            /* [out] */ double __RPC_FAR *pThirdDerivs) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamAtLength( 
            /* [in] */ double FromParam,
            /* [in] */ double Length,
            /* [out] */ double __RPC_FAR *pParam) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetLengthAtParam( 
            /* [in] */ double FromParam,
            /* [in] */ double ToParam,
            /* [out] */ double __RPC_FAR *pLength) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetContinuity( 
            /* [out] */ DWORD __RPC_FAR *nLevel) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamAnomaly( 
            /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsSingular) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetGeometryForm( 
            /* [out] */ DWORD __RPC_FAR *pForm) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMCurveVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCurve __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCurve __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCurveType )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ IID __RPC_FAR *pIID);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEndPoints )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pStartPoint[ 3 ],
            /* [out] */ double __RPC_FAR pEndPoint[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamExtents )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *pMinParam,
            /* [out] */ double __RPC_FAR *pMaxParam);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAtPoint )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPointAtParam )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTangent )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pTangents);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCurvature )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pDirections,
            /* [out] */ double __RPC_FAR *pCurvatures);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetDerivatives )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pFirstDerivs,
            /* [out] */ double __RPC_FAR *pSecondDerivs,
            /* [out] */ double __RPC_FAR *pThirdDerivs);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAtLength )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ double FromParam,
            /* [in] */ double Length,
            /* [out] */ double __RPC_FAR *pParam);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetLengthAtParam )( 
            IDMCurve __RPC_FAR * This,
            /* [in] */ double FromParam,
            /* [in] */ double ToParam,
            /* [out] */ double __RPC_FAR *pLength);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetContinuity )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *nLevel);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAnomaly )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsSingular);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetGeometryForm )( 
            IDMCurve __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *pForm);
        
        END_INTERFACE
    } IDMCurveVtbl;

    interface IDMCurve
    {
        CONST_VTBL struct IDMCurveVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCurve_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCurve_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCurve_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCurve_GetCurveType(This,pIID)	\
    (This)->lpVtbl -> GetCurveType(This,pIID)

#define IDMCurve_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#define IDMCurve_GetEndPoints(This,pStartPoint,pEndPoint)	\
    (This)->lpVtbl -> GetEndPoints(This,pStartPoint,pEndPoint)

#define IDMCurve_GetParamExtents(This,pMinParam,pMaxParam)	\
    (This)->lpVtbl -> GetParamExtents(This,pMinParam,pMaxParam)

#define IDMCurve_GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)	\
    (This)->lpVtbl -> GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)

#define IDMCurve_GetPointAtParam(This,nParams,pParams,pPoints)	\
    (This)->lpVtbl -> GetPointAtParam(This,nParams,pParams,pPoints)

#define IDMCurve_GetTangent(This,nParams,pParams,pTangents)	\
    (This)->lpVtbl -> GetTangent(This,nParams,pParams,pTangents)

#define IDMCurve_GetCurvature(This,nParams,pParams,pDirections,pCurvatures)	\
    (This)->lpVtbl -> GetCurvature(This,nParams,pParams,pDirections,pCurvatures)

#define IDMCurve_GetDerivatives(This,nParams,pParams,pFirstDerivs,pSecondDerivs,pThirdDerivs)	\
    (This)->lpVtbl -> GetDerivatives(This,nParams,pParams,pFirstDerivs,pSecondDerivs,pThirdDerivs)

#define IDMCurve_GetParamAtLength(This,FromParam,Length,pParam)	\
    (This)->lpVtbl -> GetParamAtLength(This,FromParam,Length,pParam)

#define IDMCurve_GetLengthAtParam(This,FromParam,ToParam,pLength)	\
    (This)->lpVtbl -> GetLengthAtParam(This,FromParam,ToParam,pLength)

#define IDMCurve_GetContinuity(This,nLevel)	\
    (This)->lpVtbl -> GetContinuity(This,nLevel)

#define IDMCurve_GetParamAnomaly(This,pPeriodicity,pIsSingular)	\
    (This)->lpVtbl -> GetParamAnomaly(This,pPeriodicity,pIsSingular)

#define IDMCurve_GetGeometryForm(This,pForm)	\
    (This)->lpVtbl -> GetGeometryForm(This,pForm)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCurve_GetCurveType_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ IID __RPC_FAR *pIID);


void __RPC_STUB IDMCurve_GetCurveType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetRangeBox_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 3 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 3 ]);


void __RPC_STUB IDMCurve_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetEndPoints_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pStartPoint[ 3 ],
    /* [out] */ double __RPC_FAR pEndPoint[ 3 ]);


void __RPC_STUB IDMCurve_GetEndPoints_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetParamExtents_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *pMinParam,
    /* [out] */ double __RPC_FAR *pMaxParam);


void __RPC_STUB IDMCurve_GetParamExtents_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_RemoteGetParamAtPoint_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);


void __RPC_STUB IDMCurve_RemoteGetParamAtPoint_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetPointAtParam_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in][size_is] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pPoints);


void __RPC_STUB IDMCurve_GetPointAtParam_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetTangent_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pTangents);


void __RPC_STUB IDMCurve_GetTangent_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_RemoteGetCurvature_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nDirections,
    /* [size_is][out] */ double __RPC_FAR *pDirections,
    /* [in] */ ULONG nCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pCurvatures);


void __RPC_STUB IDMCurve_RemoteGetCurvature_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_RemoteGetDerivatives_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFirstDerivs,
    /* [size_is][out] */ double __RPC_FAR *pFirstDerivs,
    /* [in] */ ULONG nSecondDerivs,
    /* [size_is][out] */ double __RPC_FAR *pSecondDerivs,
    /* [in] */ ULONG nThirdDerivs,
    /* [size_is][out] */ double __RPC_FAR *pThirdDerivs);


void __RPC_STUB IDMCurve_RemoteGetDerivatives_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetParamAtLength_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ double FromParam,
    /* [in] */ double Length,
    /* [out] */ double __RPC_FAR *pParam);


void __RPC_STUB IDMCurve_GetParamAtLength_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetLengthAtParam_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ double FromParam,
    /* [in] */ double ToParam,
    /* [out] */ double __RPC_FAR *pLength);


void __RPC_STUB IDMCurve_GetLengthAtParam_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetContinuity_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *nLevel);


void __RPC_STUB IDMCurve_GetContinuity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetParamAnomaly_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
    /* [out] */ boolean __RPC_FAR *pIsSingular);


void __RPC_STUB IDMCurve_GetParamAnomaly_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve_GetGeometryForm_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *pForm);


void __RPC_STUB IDMCurve_GetGeometryForm_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCurve_INTERFACE_DEFINED__ */


#ifndef __IDMCurve2D_INTERFACE_DEFINED__
#define __IDMCurve2D_INTERFACE_DEFINED__

/* interface IDMCurve2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCurve2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D212-0000-0000-C000-000000000046")
    IDMCurve2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetCurveType( 
            /* [out] */ IID __RPC_FAR *pIID) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetRangeBox( 
            /* [out] */ double __RPC_FAR pMinPoint[ 2 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 2 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetEndPoints( 
            /* [out] */ double __RPC_FAR pStartPoint[ 2 ],
            /* [out] */ double __RPC_FAR pEndPoint[ 2 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamExtents( 
            /* [out] */ double __RPC_FAR *pMinParam,
            /* [out] */ double __RPC_FAR *pMaxParam) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetParamAtPoint( 
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPointAtParam( 
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetTangent( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pTangents) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetCurvature( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pDirections,
            /* [out] */ double __RPC_FAR *pCurvatures) = 0;
        
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE GetDerivatives( 
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pFirstDerivs,
            /* [out] */ double __RPC_FAR *pSecondDerivs,
            /* [out] */ double __RPC_FAR *pThirdDerivs) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamAtLength( 
            /* [in] */ double FromParam,
            /* [in] */ double Length,
            /* [out] */ double __RPC_FAR *pParam) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetLengthAtParam( 
            /* [in] */ double FromParam,
            /* [in] */ double ToParam,
            /* [out] */ double __RPC_FAR *pLength) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetContinuity( 
            /* [out] */ DWORD __RPC_FAR *nLevel) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetParamAnomaly( 
            /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsSingular) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetGeometryForm( 
            /* [out] */ DWORD __RPC_FAR *pForm) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMCurve2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCurve2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCurve2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCurveType )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ IID __RPC_FAR *pIID);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetRangeBox )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pMinPoint[ 2 ],
            /* [out] */ double __RPC_FAR pMaxPoint[ 2 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEndPoints )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pStartPoint[ 2 ],
            /* [out] */ double __RPC_FAR pEndPoint[ 2 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamExtents )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *pMinParam,
            /* [out] */ double __RPC_FAR *pMaxParam);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAtPoint )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nPoints,
            /* [in] */ double __RPC_FAR *pPoints,
            /* [in] */ double __RPC_FAR *pGuessParams,
            /* [out] */ double __RPC_FAR *pMaxDeviations,
            /* [out] */ double __RPC_FAR *pParams,
            /* [out] */ DWORD __RPC_FAR *pFlags);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPointAtParam )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in][size_is] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pPoints);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTangent )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out][size_is] */ double __RPC_FAR *pTangents);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCurvature )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pDirections,
            /* [out] */ double __RPC_FAR *pCurvatures);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetDerivatives )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nParams,
            /* [in] */ double __RPC_FAR *pParams,
            /* [out] */ double __RPC_FAR *pFirstDerivs,
            /* [out] */ double __RPC_FAR *pSecondDerivs,
            /* [out] */ double __RPC_FAR *pThirdDerivs);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAtLength )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ double FromParam,
            /* [in] */ double Length,
            /* [out] */ double __RPC_FAR *pParam);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetLengthAtParam )( 
            IDMCurve2D __RPC_FAR * This,
            /* [in] */ double FromParam,
            /* [in] */ double ToParam,
            /* [out] */ double __RPC_FAR *pLength);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetContinuity )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *nLevel);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetParamAnomaly )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsSingular);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetGeometryForm )( 
            IDMCurve2D __RPC_FAR * This,
            /* [out] */ DWORD __RPC_FAR *pForm);
        
        END_INTERFACE
    } IDMCurve2DVtbl;

    interface IDMCurve2D
    {
        CONST_VTBL struct IDMCurve2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCurve2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCurve2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCurve2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCurve2D_GetCurveType(This,pIID)	\
    (This)->lpVtbl -> GetCurveType(This,pIID)

#define IDMCurve2D_GetRangeBox(This,pMinPoint,pMaxPoint)	\
    (This)->lpVtbl -> GetRangeBox(This,pMinPoint,pMaxPoint)

#define IDMCurve2D_GetEndPoints(This,pStartPoint,pEndPoint)	\
    (This)->lpVtbl -> GetEndPoints(This,pStartPoint,pEndPoint)

#define IDMCurve2D_GetParamExtents(This,pMinParam,pMaxParam)	\
    (This)->lpVtbl -> GetParamExtents(This,pMinParam,pMaxParam)

#define IDMCurve2D_GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)	\
    (This)->lpVtbl -> GetParamAtPoint(This,nPoints,pPoints,pGuessParams,pMaxDeviations,pParams,pFlags)

#define IDMCurve2D_GetPointAtParam(This,nParams,pParams,pPoints)	\
    (This)->lpVtbl -> GetPointAtParam(This,nParams,pParams,pPoints)

#define IDMCurve2D_GetTangent(This,nParams,pParams,pTangents)	\
    (This)->lpVtbl -> GetTangent(This,nParams,pParams,pTangents)

#define IDMCurve2D_GetCurvature(This,nParams,pParams,pDirections,pCurvatures)	\
    (This)->lpVtbl -> GetCurvature(This,nParams,pParams,pDirections,pCurvatures)

#define IDMCurve2D_GetDerivatives(This,nParams,pParams,pFirstDerivs,pSecondDerivs,pThirdDerivs)	\
    (This)->lpVtbl -> GetDerivatives(This,nParams,pParams,pFirstDerivs,pSecondDerivs,pThirdDerivs)

#define IDMCurve2D_GetParamAtLength(This,FromParam,Length,pParam)	\
    (This)->lpVtbl -> GetParamAtLength(This,FromParam,Length,pParam)

#define IDMCurve2D_GetLengthAtParam(This,FromParam,ToParam,pLength)	\
    (This)->lpVtbl -> GetLengthAtParam(This,FromParam,ToParam,pLength)

#define IDMCurve2D_GetContinuity(This,nLevel)	\
    (This)->lpVtbl -> GetContinuity(This,nLevel)

#define IDMCurve2D_GetParamAnomaly(This,pPeriodicity,pIsSingular)	\
    (This)->lpVtbl -> GetParamAnomaly(This,pPeriodicity,pIsSingular)

#define IDMCurve2D_GetGeometryForm(This,pForm)	\
    (This)->lpVtbl -> GetGeometryForm(This,pForm)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCurve2D_GetCurveType_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ IID __RPC_FAR *pIID);


void __RPC_STUB IDMCurve2D_GetCurveType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetRangeBox_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pMinPoint[ 2 ],
    /* [out] */ double __RPC_FAR pMaxPoint[ 2 ]);


void __RPC_STUB IDMCurve2D_GetRangeBox_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetEndPoints_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pStartPoint[ 2 ],
    /* [out] */ double __RPC_FAR pEndPoint[ 2 ]);


void __RPC_STUB IDMCurve2D_GetEndPoints_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetParamExtents_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *pMinParam,
    /* [out] */ double __RPC_FAR *pMaxParam);


void __RPC_STUB IDMCurve2D_GetParamExtents_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_RemoteGetParamAtPoint_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);


void __RPC_STUB IDMCurve2D_RemoteGetParamAtPoint_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetPointAtParam_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in][size_is] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pPoints);


void __RPC_STUB IDMCurve2D_GetPointAtParam_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetTangent_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out][size_is] */ double __RPC_FAR *pTangents);


void __RPC_STUB IDMCurve2D_GetTangent_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_RemoteGetCurvature_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nDirections,
    /* [size_is][out] */ double __RPC_FAR *pDirections,
    /* [in] */ ULONG nCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pCurvatures);


void __RPC_STUB IDMCurve2D_RemoteGetCurvature_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_RemoteGetDerivatives_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFirstDerivs,
    /* [size_is][out] */ double __RPC_FAR *pFirstDerivs,
    /* [in] */ ULONG nSecondDerivs,
    /* [size_is][out] */ double __RPC_FAR *pSecondDerivs,
    /* [in] */ ULONG nThirdDerivs,
    /* [size_is][out] */ double __RPC_FAR *pThirdDerivs);


void __RPC_STUB IDMCurve2D_RemoteGetDerivatives_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetParamAtLength_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ double FromParam,
    /* [in] */ double Length,
    /* [out] */ double __RPC_FAR *pParam);


void __RPC_STUB IDMCurve2D_GetParamAtLength_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetLengthAtParam_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ double FromParam,
    /* [in] */ double ToParam,
    /* [out] */ double __RPC_FAR *pLength);


void __RPC_STUB IDMCurve2D_GetLengthAtParam_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetContinuity_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *nLevel);


void __RPC_STUB IDMCurve2D_GetContinuity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetParamAnomaly_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pPeriodicity[ 2 ],
    /* [out] */ boolean __RPC_FAR *pIsSingular);


void __RPC_STUB IDMCurve2D_GetParamAnomaly_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMCurve2D_GetGeometryForm_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [out] */ DWORD __RPC_FAR *pForm);


void __RPC_STUB IDMCurve2D_GetGeometryForm_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCurve2D_INTERFACE_DEFINED__ */


#ifndef __IDMCone_INTERFACE_DEFINED__
#define __IDMCone_INTERFACE_DEFINED__

/* interface IDMCone */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCone;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D213-0000-0000-C000-000000000046")
    IDMCone : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetConeData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius,
            /* [out] */ double __RPC_FAR *pHalfAngle,
            /* [out] */ boolean __RPC_FAR *pIsExpanding) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMConeVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCone __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCone __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCone __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetConeData )( 
            IDMCone __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius,
            /* [out] */ double __RPC_FAR *pHalfAngle,
            /* [out] */ boolean __RPC_FAR *pIsExpanding);
        
        END_INTERFACE
    } IDMConeVtbl;

    interface IDMCone
    {
        CONST_VTBL struct IDMConeVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCone_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCone_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCone_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCone_GetConeData(This,pBasePoint,pAxisVector,pRadius,pHalfAngle,pIsExpanding)	\
    (This)->lpVtbl -> GetConeData(This,pBasePoint,pAxisVector,pRadius,pHalfAngle,pIsExpanding)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCone_GetConeData_Proxy( 
    IDMCone __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
    /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
    /* [out] */ double __RPC_FAR *pRadius,
    /* [out] */ double __RPC_FAR *pHalfAngle,
    /* [out] */ boolean __RPC_FAR *pIsExpanding);


void __RPC_STUB IDMCone_GetConeData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCone_INTERFACE_DEFINED__ */


#ifndef __IDMCylinder_INTERFACE_DEFINED__
#define __IDMCylinder_INTERFACE_DEFINED__

/* interface IDMCylinder */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCylinder;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D214-0000-0000-C000-000000000046")
    IDMCylinder : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetCylinderData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMCylinderVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCylinder __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCylinder __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCylinder __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCylinderData )( 
            IDMCylinder __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius);
        
        END_INTERFACE
    } IDMCylinderVtbl;

    interface IDMCylinder
    {
        CONST_VTBL struct IDMCylinderVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCylinder_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCylinder_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCylinder_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCylinder_GetCylinderData(This,pBasePoint,pAxisVector,pRadius)	\
    (This)->lpVtbl -> GetCylinderData(This,pBasePoint,pAxisVector,pRadius)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCylinder_GetCylinderData_Proxy( 
    IDMCylinder __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
    /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
    /* [out] */ double __RPC_FAR *pRadius);


void __RPC_STUB IDMCylinder_GetCylinderData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCylinder_INTERFACE_DEFINED__ */


#ifndef __IDMSphere_INTERFACE_DEFINED__
#define __IDMSphere_INTERFACE_DEFINED__

/* interface IDMSphere */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMSphere;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D215-0000-0000-C000-000000000046")
    IDMSphere : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetSphereData( 
            /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMSphereVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMSphere __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMSphere __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMSphere __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetSphereData )( 
            IDMSphere __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius);
        
        END_INTERFACE
    } IDMSphereVtbl;

    interface IDMSphere
    {
        CONST_VTBL struct IDMSphereVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMSphere_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMSphere_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMSphere_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMSphere_GetSphereData(This,pCenterPoint,pRadius)	\
    (This)->lpVtbl -> GetSphereData(This,pCenterPoint,pRadius)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMSphere_GetSphereData_Proxy( 
    IDMSphere __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
    /* [out] */ double __RPC_FAR *pRadius);


void __RPC_STUB IDMSphere_GetSphereData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMSphere_INTERFACE_DEFINED__ */


#ifndef __IDMTorus_INTERFACE_DEFINED__
#define __IDMTorus_INTERFACE_DEFINED__

/* interface IDMTorus */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMTorus;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D216-0000-0000-C000-000000000046")
    IDMTorus : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetTorusData( 
            /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pMajorRadius,
            /* [out] */ double __RPC_FAR *pMinorRadius) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMTorusVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMTorus __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMTorus __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMTorus __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTorusData )( 
            IDMTorus __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pMajorRadius,
            /* [out] */ double __RPC_FAR *pMinorRadius);
        
        END_INTERFACE
    } IDMTorusVtbl;

    interface IDMTorus
    {
        CONST_VTBL struct IDMTorusVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMTorus_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMTorus_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMTorus_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMTorus_GetTorusData(This,pCenterPoint,pAxisVector,pMajorRadius,pMinorRadius)	\
    (This)->lpVtbl -> GetTorusData(This,pCenterPoint,pAxisVector,pMajorRadius,pMinorRadius)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMTorus_GetTorusData_Proxy( 
    IDMTorus __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pCenterPoint[ 3 ],
    /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
    /* [out] */ double __RPC_FAR *pMajorRadius,
    /* [out] */ double __RPC_FAR *pMinorRadius);


void __RPC_STUB IDMTorus_GetTorusData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMTorus_INTERFACE_DEFINED__ */


#ifndef __IDMPlane_INTERFACE_DEFINED__
#define __IDMPlane_INTERFACE_DEFINED__

/* interface IDMPlane */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMPlane;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D217-0000-0000-C000-000000000046")
    IDMPlane : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetPlaneData( 
            /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
            /* [out] */ double __RPC_FAR pNormalVector[ 3 ]) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMPlaneVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMPlane __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMPlane __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMPlane __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPlaneData )( 
            IDMPlane __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
            /* [out] */ double __RPC_FAR pNormalVector[ 3 ]);
        
        END_INTERFACE
    } IDMPlaneVtbl;

    interface IDMPlane
    {
        CONST_VTBL struct IDMPlaneVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMPlane_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMPlane_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMPlane_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMPlane_GetPlaneData(This,pRootPoint,pNormalVector)	\
    (This)->lpVtbl -> GetPlaneData(This,pRootPoint,pNormalVector)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMPlane_GetPlaneData_Proxy( 
    IDMPlane __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
    /* [out] */ double __RPC_FAR pNormalVector[ 3 ]);


void __RPC_STUB IDMPlane_GetPlaneData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMPlane_INTERFACE_DEFINED__ */


#ifndef __IDMBSplineSurface_INTERFACE_DEFINED__
#define __IDMBSplineSurface_INTERFACE_DEFINED__

/* interface IDMBSplineSurface */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMBSplineSurface;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D218-0000-0000-C000-000000000046")
    IDMBSplineSurface : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetBSplineInfo( 
            /* [out] */ ULONG __RPC_FAR pOrder[ 2 ],
            /* [out] */ ULONG __RPC_FAR pNumPoles[ 2 ],
            /* [out] */ ULONG __RPC_FAR pNumKnots[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR pIsClosed[ 2 ],
            /* [out] */ boolean __RPC_FAR pIsPeriodic[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsPlanar) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetBSplineData( 
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnotsU,
            /* [in] */ ULONG nKnotsV,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnotsU,
            /* [out][size_is] */ double __RPC_FAR *pKnotsV,
            /* [out][size_is] */ double __RPC_FAR *pWeights) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMBSplineSurfaceVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMBSplineSurface __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMBSplineSurface __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMBSplineSurface __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineInfo )( 
            IDMBSplineSurface __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR pOrder[ 2 ],
            /* [out] */ ULONG __RPC_FAR pNumPoles[ 2 ],
            /* [out] */ ULONG __RPC_FAR pNumKnots[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR pIsClosed[ 2 ],
            /* [out] */ boolean __RPC_FAR pIsPeriodic[ 2 ],
            /* [out] */ boolean __RPC_FAR *pIsPlanar);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineData )( 
            IDMBSplineSurface __RPC_FAR * This,
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnotsU,
            /* [in] */ ULONG nKnotsV,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnotsU,
            /* [out][size_is] */ double __RPC_FAR *pKnotsV,
            /* [out][size_is] */ double __RPC_FAR *pWeights);
        
        END_INTERFACE
    } IDMBSplineSurfaceVtbl;

    interface IDMBSplineSurface
    {
        CONST_VTBL struct IDMBSplineSurfaceVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMBSplineSurface_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMBSplineSurface_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMBSplineSurface_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMBSplineSurface_GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic,pIsPlanar)	\
    (This)->lpVtbl -> GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic,pIsPlanar)

#define IDMBSplineSurface_GetBSplineData(This,nPoles,nKnotsU,nKnotsV,nWeights,pPoles,pKnotsU,pKnotsV,pWeights)	\
    (This)->lpVtbl -> GetBSplineData(This,nPoles,nKnotsU,nKnotsV,nWeights,pPoles,pKnotsU,pKnotsV,pWeights)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMBSplineSurface_GetBSplineInfo_Proxy( 
    IDMBSplineSurface __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR pOrder[ 2 ],
    /* [out] */ ULONG __RPC_FAR pNumPoles[ 2 ],
    /* [out] */ ULONG __RPC_FAR pNumKnots[ 2 ],
    /* [out] */ boolean __RPC_FAR *pIsRational,
    /* [out] */ boolean __RPC_FAR pIsClosed[ 2 ],
    /* [out] */ boolean __RPC_FAR pIsPeriodic[ 2 ],
    /* [out] */ boolean __RPC_FAR *pIsPlanar);


void __RPC_STUB IDMBSplineSurface_GetBSplineInfo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMBSplineSurface_GetBSplineData_Proxy( 
    IDMBSplineSurface __RPC_FAR * This,
    /* [in] */ ULONG nPoles,
    /* [in] */ ULONG nKnotsU,
    /* [in] */ ULONG nKnotsV,
    /* [in] */ ULONG nWeights,
    /* [out][size_is] */ double __RPC_FAR *pPoles,
    /* [out][size_is] */ double __RPC_FAR *pKnotsU,
    /* [out][size_is] */ double __RPC_FAR *pKnotsV,
    /* [out][size_is] */ double __RPC_FAR *pWeights);


void __RPC_STUB IDMBSplineSurface_GetBSplineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMBSplineSurface_INTERFACE_DEFINED__ */


#ifndef __IDMCircle_INTERFACE_DEFINED__
#define __IDMCircle_INTERFACE_DEFINED__

/* interface IDMCircle */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCircle;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D219-0000-0000-C000-000000000046")
    IDMCircle : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetCircleData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMCircleVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCircle __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCircle __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCircle __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCircleData )( 
            IDMCircle __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR *pRadius);
        
        END_INTERFACE
    } IDMCircleVtbl;

    interface IDMCircle
    {
        CONST_VTBL struct IDMCircleVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCircle_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCircle_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCircle_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCircle_GetCircleData(This,pBasePoint,pAxisVector,pRadius)	\
    (This)->lpVtbl -> GetCircleData(This,pBasePoint,pAxisVector,pRadius)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCircle_GetCircleData_Proxy( 
    IDMCircle __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
    /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
    /* [out] */ double __RPC_FAR *pRadius);


void __RPC_STUB IDMCircle_GetCircleData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCircle_INTERFACE_DEFINED__ */


#ifndef __IDMEllipse_INTERFACE_DEFINED__
#define __IDMEllipse_INTERFACE_DEFINED__

/* interface IDMEllipse */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMEllipse;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21A-0000-0000-C000-000000000046")
    IDMEllipse : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetEllipseData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR pMajorAxis[ 3 ],
            /* [out] */ double __RPC_FAR *pMinorMajorRatio) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMEllipseVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMEllipse __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMEllipse __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMEllipse __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEllipseData )( 
            IDMEllipse __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
            /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
            /* [out] */ double __RPC_FAR pMajorAxis[ 3 ],
            /* [out] */ double __RPC_FAR *pMinorMajorRatio);
        
        END_INTERFACE
    } IDMEllipseVtbl;

    interface IDMEllipse
    {
        CONST_VTBL struct IDMEllipseVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMEllipse_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMEllipse_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMEllipse_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMEllipse_GetEllipseData(This,pBasePoint,pAxisVector,pMajorAxis,pMinorMajorRatio)	\
    (This)->lpVtbl -> GetEllipseData(This,pBasePoint,pAxisVector,pMajorAxis,pMinorMajorRatio)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMEllipse_GetEllipseData_Proxy( 
    IDMEllipse __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 3 ],
    /* [out] */ double __RPC_FAR pAxisVector[ 3 ],
    /* [out] */ double __RPC_FAR pMajorAxis[ 3 ],
    /* [out] */ double __RPC_FAR *pMinorMajorRatio);


void __RPC_STUB IDMEllipse_GetEllipseData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMEllipse_INTERFACE_DEFINED__ */


#ifndef __IDMLine_INTERFACE_DEFINED__
#define __IDMLine_INTERFACE_DEFINED__

/* interface IDMLine */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMLine;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21B-0000-0000-C000-000000000046")
    IDMLine : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetLineData( 
            /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
            /* [out] */ double __RPC_FAR pDirectionVector[ 3 ]) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMLineVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMLine __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMLine __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMLine __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetLineData )( 
            IDMLine __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
            /* [out] */ double __RPC_FAR pDirectionVector[ 3 ]);
        
        END_INTERFACE
    } IDMLineVtbl;

    interface IDMLine
    {
        CONST_VTBL struct IDMLineVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMLine_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMLine_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMLine_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMLine_GetLineData(This,pRootPoint,pDirectionVector)	\
    (This)->lpVtbl -> GetLineData(This,pRootPoint,pDirectionVector)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMLine_GetLineData_Proxy( 
    IDMLine __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pRootPoint[ 3 ],
    /* [out] */ double __RPC_FAR pDirectionVector[ 3 ]);


void __RPC_STUB IDMLine_GetLineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMLine_INTERFACE_DEFINED__ */


#ifndef __IDMPolyline_INTERFACE_DEFINED__
#define __IDMPolyline_INTERFACE_DEFINED__

/* interface IDMPolyline */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMPolyline;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21C-0000-0000-C000-000000000046")
    IDMPolyline : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetPolylineInfo( 
            /* [out] */ ULONG __RPC_FAR *pNumPoints,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPlanar,
            /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPolylineData( 
            /* [in] */ ULONG nPoints,
            /* [out][size_is] */ double __RPC_FAR *pPoints) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMPolylineVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMPolyline __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMPolyline __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMPolyline __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPolylineInfo )( 
            IDMPolyline __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pNumPoints,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPlanar,
            /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPolylineData )( 
            IDMPolyline __RPC_FAR * This,
            /* [in] */ ULONG nPoints,
            /* [out][size_is] */ double __RPC_FAR *pPoints);
        
        END_INTERFACE
    } IDMPolylineVtbl;

    interface IDMPolyline
    {
        CONST_VTBL struct IDMPolylineVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMPolyline_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMPolyline_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMPolyline_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMPolyline_GetPolylineInfo(This,pNumPoints,pIsClosed,pIsPlanar,pPlaneVector)	\
    (This)->lpVtbl -> GetPolylineInfo(This,pNumPoints,pIsClosed,pIsPlanar,pPlaneVector)

#define IDMPolyline_GetPolylineData(This,nPoints,pPoints)	\
    (This)->lpVtbl -> GetPolylineData(This,nPoints,pPoints)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMPolyline_GetPolylineInfo_Proxy( 
    IDMPolyline __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pNumPoints,
    /* [out] */ boolean __RPC_FAR *pIsClosed,
    /* [out] */ boolean __RPC_FAR *pIsPlanar,
    /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]);


void __RPC_STUB IDMPolyline_GetPolylineInfo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMPolyline_GetPolylineData_Proxy( 
    IDMPolyline __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [out][size_is] */ double __RPC_FAR *pPoints);


void __RPC_STUB IDMPolyline_GetPolylineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMPolyline_INTERFACE_DEFINED__ */


#ifndef __IDMBSplineCurve_INTERFACE_DEFINED__
#define __IDMBSplineCurve_INTERFACE_DEFINED__

/* interface IDMBSplineCurve */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMBSplineCurve;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21D-0000-0000-C000-000000000046")
    IDMBSplineCurve : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetBSplineInfo( 
            /* [out] */ ULONG __RPC_FAR *pOrder,
            /* [out] */ ULONG __RPC_FAR *pNumPoles,
            /* [out] */ ULONG __RPC_FAR *pNumKnots,
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPeriodic,
            /* [out] */ boolean __RPC_FAR *pIsPlanar,
            /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetBSplineData( 
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnots,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnots,
            /* [out][size_is] */ double __RPC_FAR *pWeights) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMBSplineCurveVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMBSplineCurve __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMBSplineCurve __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMBSplineCurve __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineInfo )( 
            IDMBSplineCurve __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pOrder,
            /* [out] */ ULONG __RPC_FAR *pNumPoles,
            /* [out] */ ULONG __RPC_FAR *pNumKnots,
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPeriodic,
            /* [out] */ boolean __RPC_FAR *pIsPlanar,
            /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineData )( 
            IDMBSplineCurve __RPC_FAR * This,
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnots,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnots,
            /* [out][size_is] */ double __RPC_FAR *pWeights);
        
        END_INTERFACE
    } IDMBSplineCurveVtbl;

    interface IDMBSplineCurve
    {
        CONST_VTBL struct IDMBSplineCurveVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMBSplineCurve_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMBSplineCurve_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMBSplineCurve_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMBSplineCurve_GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic,pIsPlanar,pPlaneVector)	\
    (This)->lpVtbl -> GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic,pIsPlanar,pPlaneVector)

#define IDMBSplineCurve_GetBSplineData(This,nPoles,nKnots,nWeights,pPoles,pKnots,pWeights)	\
    (This)->lpVtbl -> GetBSplineData(This,nPoles,nKnots,nWeights,pPoles,pKnots,pWeights)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMBSplineCurve_GetBSplineInfo_Proxy( 
    IDMBSplineCurve __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pOrder,
    /* [out] */ ULONG __RPC_FAR *pNumPoles,
    /* [out] */ ULONG __RPC_FAR *pNumKnots,
    /* [out] */ boolean __RPC_FAR *pIsRational,
    /* [out] */ boolean __RPC_FAR *pIsClosed,
    /* [out] */ boolean __RPC_FAR *pIsPeriodic,
    /* [out] */ boolean __RPC_FAR *pIsPlanar,
    /* [out] */ double __RPC_FAR pPlaneVector[ 3 ]);


void __RPC_STUB IDMBSplineCurve_GetBSplineInfo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMBSplineCurve_GetBSplineData_Proxy( 
    IDMBSplineCurve __RPC_FAR * This,
    /* [in] */ ULONG nPoles,
    /* [in] */ ULONG nKnots,
    /* [in] */ ULONG nWeights,
    /* [out][size_is] */ double __RPC_FAR *pPoles,
    /* [out][size_is] */ double __RPC_FAR *pKnots,
    /* [out][size_is] */ double __RPC_FAR *pWeights);


void __RPC_STUB IDMBSplineCurve_GetBSplineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMBSplineCurve_INTERFACE_DEFINED__ */


#ifndef __IDMCircle2D_INTERFACE_DEFINED__
#define __IDMCircle2D_INTERFACE_DEFINED__

/* interface IDMCircle2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMCircle2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21E-0000-0000-C000-000000000046")
    IDMCircle2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetCircleData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
            /* [out] */ double __RPC_FAR *pRadius) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMCircle2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMCircle2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMCircle2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMCircle2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCircleData )( 
            IDMCircle2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
            /* [out] */ double __RPC_FAR *pRadius);
        
        END_INTERFACE
    } IDMCircle2DVtbl;

    interface IDMCircle2D
    {
        CONST_VTBL struct IDMCircle2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMCircle2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMCircle2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMCircle2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMCircle2D_GetCircleData(This,pBasePoint,pRadius)	\
    (This)->lpVtbl -> GetCircleData(This,pBasePoint,pRadius)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMCircle2D_GetCircleData_Proxy( 
    IDMCircle2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
    /* [out] */ double __RPC_FAR *pRadius);


void __RPC_STUB IDMCircle2D_GetCircleData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMCircle2D_INTERFACE_DEFINED__ */


#ifndef __IDMEllipse2D_INTERFACE_DEFINED__
#define __IDMEllipse2D_INTERFACE_DEFINED__

/* interface IDMEllipse2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMEllipse2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D21F-0000-0000-C000-000000000046")
    IDMEllipse2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetEllipseData( 
            /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
            /* [out] */ double __RPC_FAR pMajorAxis[ 2 ],
            /* [out] */ double __RPC_FAR *pMinorMajorRatio) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMEllipse2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMEllipse2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMEllipse2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMEllipse2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetEllipseData )( 
            IDMEllipse2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
            /* [out] */ double __RPC_FAR pMajorAxis[ 2 ],
            /* [out] */ double __RPC_FAR *pMinorMajorRatio);
        
        END_INTERFACE
    } IDMEllipse2DVtbl;

    interface IDMEllipse2D
    {
        CONST_VTBL struct IDMEllipse2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMEllipse2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMEllipse2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMEllipse2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMEllipse2D_GetEllipseData(This,pBasePoint,pMajorAxis,pMinorMajorRatio)	\
    (This)->lpVtbl -> GetEllipseData(This,pBasePoint,pMajorAxis,pMinorMajorRatio)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMEllipse2D_GetEllipseData_Proxy( 
    IDMEllipse2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pBasePoint[ 2 ],
    /* [out] */ double __RPC_FAR pMajorAxis[ 2 ],
    /* [out] */ double __RPC_FAR *pMinorMajorRatio);


void __RPC_STUB IDMEllipse2D_GetEllipseData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMEllipse2D_INTERFACE_DEFINED__ */


#ifndef __IDMLine2D_INTERFACE_DEFINED__
#define __IDMLine2D_INTERFACE_DEFINED__

/* interface IDMLine2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMLine2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D220-0000-0000-C000-000000000046")
    IDMLine2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetLineData( 
            /* [out] */ double __RPC_FAR pRootPoint[ 2 ],
            /* [out] */ double __RPC_FAR pDirectionVector[ 2 ]) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMLine2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMLine2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMLine2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMLine2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetLineData )( 
            IDMLine2D __RPC_FAR * This,
            /* [out] */ double __RPC_FAR pRootPoint[ 2 ],
            /* [out] */ double __RPC_FAR pDirectionVector[ 2 ]);
        
        END_INTERFACE
    } IDMLine2DVtbl;

    interface IDMLine2D
    {
        CONST_VTBL struct IDMLine2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMLine2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMLine2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMLine2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMLine2D_GetLineData(This,pRootPoint,pDirectionVector)	\
    (This)->lpVtbl -> GetLineData(This,pRootPoint,pDirectionVector)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMLine2D_GetLineData_Proxy( 
    IDMLine2D __RPC_FAR * This,
    /* [out] */ double __RPC_FAR pRootPoint[ 2 ],
    /* [out] */ double __RPC_FAR pDirectionVector[ 2 ]);


void __RPC_STUB IDMLine2D_GetLineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMLine2D_INTERFACE_DEFINED__ */


#ifndef __IDMPolyline2D_INTERFACE_DEFINED__
#define __IDMPolyline2D_INTERFACE_DEFINED__

/* interface IDMPolyline2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMPolyline2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D221-0000-0000-C000-000000000046")
    IDMPolyline2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetPolylineInfo( 
            /* [out] */ ULONG __RPC_FAR *pNumPoints,
            /* [out] */ boolean __RPC_FAR *pIsClosed) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPolylineData( 
            /* [in] */ ULONG nPoints,
            /* [out][size_is] */ double __RPC_FAR *pPoints) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMPolyline2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMPolyline2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMPolyline2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMPolyline2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPolylineInfo )( 
            IDMPolyline2D __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pNumPoints,
            /* [out] */ boolean __RPC_FAR *pIsClosed);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetPolylineData )( 
            IDMPolyline2D __RPC_FAR * This,
            /* [in] */ ULONG nPoints,
            /* [out][size_is] */ double __RPC_FAR *pPoints);
        
        END_INTERFACE
    } IDMPolyline2DVtbl;

    interface IDMPolyline2D
    {
        CONST_VTBL struct IDMPolyline2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMPolyline2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMPolyline2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMPolyline2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMPolyline2D_GetPolylineInfo(This,pNumPoints,pIsClosed)	\
    (This)->lpVtbl -> GetPolylineInfo(This,pNumPoints,pIsClosed)

#define IDMPolyline2D_GetPolylineData(This,nPoints,pPoints)	\
    (This)->lpVtbl -> GetPolylineData(This,nPoints,pPoints)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMPolyline2D_GetPolylineInfo_Proxy( 
    IDMPolyline2D __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pNumPoints,
    /* [out] */ boolean __RPC_FAR *pIsClosed);


void __RPC_STUB IDMPolyline2D_GetPolylineInfo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMPolyline2D_GetPolylineData_Proxy( 
    IDMPolyline2D __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [out][size_is] */ double __RPC_FAR *pPoints);


void __RPC_STUB IDMPolyline2D_GetPolylineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMPolyline2D_INTERFACE_DEFINED__ */


#ifndef __IDMBSplineCurve2D_INTERFACE_DEFINED__
#define __IDMBSplineCurve2D_INTERFACE_DEFINED__

/* interface IDMBSplineCurve2D */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMBSplineCurve2D;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D222-0000-0000-C000-000000000046")
    IDMBSplineCurve2D : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetBSplineInfo( 
            /* [out] */ ULONG __RPC_FAR *pOrder,
            /* [out] */ ULONG __RPC_FAR *pNumPoles,
            /* [out] */ ULONG __RPC_FAR *pNumKnots,
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPeriodic) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetBSplineData( 
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnots,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnots,
            /* [out][size_is] */ double __RPC_FAR *pWeights) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMBSplineCurve2DVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMBSplineCurve2D __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMBSplineCurve2D __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMBSplineCurve2D __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineInfo )( 
            IDMBSplineCurve2D __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pOrder,
            /* [out] */ ULONG __RPC_FAR *pNumPoles,
            /* [out] */ ULONG __RPC_FAR *pNumKnots,
            /* [out] */ boolean __RPC_FAR *pIsRational,
            /* [out] */ boolean __RPC_FAR *pIsClosed,
            /* [out] */ boolean __RPC_FAR *pIsPeriodic);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetBSplineData )( 
            IDMBSplineCurve2D __RPC_FAR * This,
            /* [in] */ ULONG nPoles,
            /* [in] */ ULONG nKnots,
            /* [in] */ ULONG nWeights,
            /* [out][size_is] */ double __RPC_FAR *pPoles,
            /* [out][size_is] */ double __RPC_FAR *pKnots,
            /* [out][size_is] */ double __RPC_FAR *pWeights);
        
        END_INTERFACE
    } IDMBSplineCurve2DVtbl;

    interface IDMBSplineCurve2D
    {
        CONST_VTBL struct IDMBSplineCurve2DVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMBSplineCurve2D_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMBSplineCurve2D_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMBSplineCurve2D_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMBSplineCurve2D_GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic)	\
    (This)->lpVtbl -> GetBSplineInfo(This,pOrder,pNumPoles,pNumKnots,pIsRational,pIsClosed,pIsPeriodic)

#define IDMBSplineCurve2D_GetBSplineData(This,nPoles,nKnots,nWeights,pPoles,pKnots,pWeights)	\
    (This)->lpVtbl -> GetBSplineData(This,nPoles,nKnots,nWeights,pPoles,pKnots,pWeights)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMBSplineCurve2D_GetBSplineInfo_Proxy( 
    IDMBSplineCurve2D __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pOrder,
    /* [out] */ ULONG __RPC_FAR *pNumPoles,
    /* [out] */ ULONG __RPC_FAR *pNumKnots,
    /* [out] */ boolean __RPC_FAR *pIsRational,
    /* [out] */ boolean __RPC_FAR *pIsClosed,
    /* [out] */ boolean __RPC_FAR *pIsPeriodic);


void __RPC_STUB IDMBSplineCurve2D_GetBSplineInfo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMBSplineCurve2D_GetBSplineData_Proxy( 
    IDMBSplineCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nPoles,
    /* [in] */ ULONG nKnots,
    /* [in] */ ULONG nWeights,
    /* [out][size_is] */ double __RPC_FAR *pPoles,
    /* [out][size_is] */ double __RPC_FAR *pKnots,
    /* [out][size_is] */ double __RPC_FAR *pWeights);


void __RPC_STUB IDMBSplineCurve2D_GetBSplineData_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMBSplineCurve2D_INTERFACE_DEFINED__ */


#ifndef __IEnumDMSurfaceBodies_INTERFACE_DEFINED__
#define __IEnumDMSurfaceBodies_INTERFACE_DEFINED__

/* interface IEnumDMSurfaceBodies */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMSurfaceBodies;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D223-0000-0000-C000-000000000046")
    IEnumDMSurfaceBodies : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cSurfaceBody,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pSurfaceBody,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cCurves) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMSurfaceBodiesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMSurfaceBodies __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMSurfaceBodies __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMSurfaceBodies __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMSurfaceBodies __RPC_FAR * This,
            /* [in] */ ULONG cSurfaceBody,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *pSurfaceBody,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMSurfaceBodies __RPC_FAR * This,
            /* [in] */ ULONG cCurves);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMSurfaceBodies __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMSurfaceBodies __RPC_FAR * This,
            /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies);
        
        END_INTERFACE
    } IEnumDMSurfaceBodiesVtbl;

    interface IEnumDMSurfaceBodies
    {
        CONST_VTBL struct IEnumDMSurfaceBodiesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMSurfaceBodies_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMSurfaceBodies_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMSurfaceBodies_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMSurfaceBodies_Next(This,cSurfaceBody,pSurfaceBody,pcFetched)	\
    (This)->lpVtbl -> Next(This,cSurfaceBody,pSurfaceBody,pcFetched)

#define IEnumDMSurfaceBodies_Skip(This,cCurves)	\
    (This)->lpVtbl -> Skip(This,cCurves)

#define IEnumDMSurfaceBodies_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMSurfaceBodies_Clone(This,pEnumSurfaceBodies)	\
    (This)->lpVtbl -> Clone(This,pEnumSurfaceBodies)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_RemoteNext_Proxy( 
    IEnumDMSurfaceBodies __RPC_FAR * This,
    /* [in] */ ULONG cSurfaceBody,
    /* [length_is][size_is][out] */ LPDMSURFACEBODY __RPC_FAR *pSurfaceBody,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMSurfaceBodies_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_Skip_Proxy( 
    IEnumDMSurfaceBodies __RPC_FAR * This,
    /* [in] */ ULONG cCurves);


void __RPC_STUB IEnumDMSurfaceBodies_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_Reset_Proxy( 
    IEnumDMSurfaceBodies __RPC_FAR * This);


void __RPC_STUB IEnumDMSurfaceBodies_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_Clone_Proxy( 
    IEnumDMSurfaceBodies __RPC_FAR * This,
    /* [out] */ LPENUM_DMSURFACEBODIES __RPC_FAR *pEnumSurfaceBodies);


void __RPC_STUB IEnumDMSurfaceBodies_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMSurfaceBodies_INTERFACE_DEFINED__ */


#ifndef __IEnumDMShells_INTERFACE_DEFINED__
#define __IEnumDMShells_INTERFACE_DEFINED__

/* interface IEnumDMShells */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMShells;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D224-0000-0000-C000-000000000046")
    IEnumDMShells : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cShells,
            /* [out] */ LPDMSHELL __RPC_FAR *pShell,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cShells) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMShellsVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMShells __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMShells __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMShells __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMShells __RPC_FAR * This,
            /* [in] */ ULONG cShells,
            /* [out] */ LPDMSHELL __RPC_FAR *pShell,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMShells __RPC_FAR * This,
            /* [in] */ ULONG cShells);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMShells __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMShells __RPC_FAR * This,
            /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells);
        
        END_INTERFACE
    } IEnumDMShellsVtbl;

    interface IEnumDMShells
    {
        CONST_VTBL struct IEnumDMShellsVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMShells_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMShells_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMShells_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMShells_Next(This,cShells,pShell,pcFetched)	\
    (This)->lpVtbl -> Next(This,cShells,pShell,pcFetched)

#define IEnumDMShells_Skip(This,cShells)	\
    (This)->lpVtbl -> Skip(This,cShells)

#define IEnumDMShells_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMShells_Clone(This,pEnumShells)	\
    (This)->lpVtbl -> Clone(This,pEnumShells)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMShells_RemoteNext_Proxy( 
    IEnumDMShells __RPC_FAR * This,
    /* [in] */ ULONG cShells,
    /* [length_is][size_is][out] */ LPDMSHELL __RPC_FAR *pShell,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMShells_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMShells_Skip_Proxy( 
    IEnumDMShells __RPC_FAR * This,
    /* [in] */ ULONG cShells);


void __RPC_STUB IEnumDMShells_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMShells_Reset_Proxy( 
    IEnumDMShells __RPC_FAR * This);


void __RPC_STUB IEnumDMShells_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMShells_Clone_Proxy( 
    IEnumDMShells __RPC_FAR * This,
    /* [out] */ LPENUM_DMSHELLS __RPC_FAR *pEnumShells);


void __RPC_STUB IEnumDMShells_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMShells_INTERFACE_DEFINED__ */


#ifndef __IEnumDMFaces_INTERFACE_DEFINED__
#define __IEnumDMFaces_INTERFACE_DEFINED__

/* interface IEnumDMFaces */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMFaces;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D225-0000-0000-C000-000000000046")
    IEnumDMFaces : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cFaces,
            /* [out] */ LPDMFACE __RPC_FAR *pFace,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cFaces) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMFacesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMFaces __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMFaces __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMFaces __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMFaces __RPC_FAR * This,
            /* [in] */ ULONG cFaces,
            /* [out] */ LPDMFACE __RPC_FAR *pFace,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMFaces __RPC_FAR * This,
            /* [in] */ ULONG cFaces);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMFaces __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMFaces __RPC_FAR * This,
            /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);
        
        END_INTERFACE
    } IEnumDMFacesVtbl;

    interface IEnumDMFaces
    {
        CONST_VTBL struct IEnumDMFacesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMFaces_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMFaces_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMFaces_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMFaces_Next(This,cFaces,pFace,pcFetched)	\
    (This)->lpVtbl -> Next(This,cFaces,pFace,pcFetched)

#define IEnumDMFaces_Skip(This,cFaces)	\
    (This)->lpVtbl -> Skip(This,cFaces)

#define IEnumDMFaces_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMFaces_Clone(This,pEnumFaces)	\
    (This)->lpVtbl -> Clone(This,pEnumFaces)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMFaces_RemoteNext_Proxy( 
    IEnumDMFaces __RPC_FAR * This,
    /* [in] */ ULONG cFaces,
    /* [length_is][size_is][out] */ LPDMFACE __RPC_FAR *pFace,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMFaces_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMFaces_Skip_Proxy( 
    IEnumDMFaces __RPC_FAR * This,
    /* [in] */ ULONG cFaces);


void __RPC_STUB IEnumDMFaces_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMFaces_Reset_Proxy( 
    IEnumDMFaces __RPC_FAR * This);


void __RPC_STUB IEnumDMFaces_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMFaces_Clone_Proxy( 
    IEnumDMFaces __RPC_FAR * This,
    /* [out] */ LPENUM_DMFACES __RPC_FAR *pEnumFaces);


void __RPC_STUB IEnumDMFaces_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMFaces_INTERFACE_DEFINED__ */


#ifndef __IEnumDMLoops_INTERFACE_DEFINED__
#define __IEnumDMLoops_INTERFACE_DEFINED__

/* interface IEnumDMLoops */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMLoops;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D226-0000-0000-C000-000000000046")
    IEnumDMLoops : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cEdge,
            /* [out] */ LPDMLOOP __RPC_FAR *pLoop,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cLoops) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMLoopsVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMLoops __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMLoops __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMLoops __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMLoops __RPC_FAR * This,
            /* [in] */ ULONG cEdge,
            /* [out] */ LPDMLOOP __RPC_FAR *pLoop,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMLoops __RPC_FAR * This,
            /* [in] */ ULONG cLoops);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMLoops __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMLoops __RPC_FAR * This,
            /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops);
        
        END_INTERFACE
    } IEnumDMLoopsVtbl;

    interface IEnumDMLoops
    {
        CONST_VTBL struct IEnumDMLoopsVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMLoops_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMLoops_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMLoops_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMLoops_Next(This,cEdge,pLoop,pcFetched)	\
    (This)->lpVtbl -> Next(This,cEdge,pLoop,pcFetched)

#define IEnumDMLoops_Skip(This,cLoops)	\
    (This)->lpVtbl -> Skip(This,cLoops)

#define IEnumDMLoops_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMLoops_Clone(This,pEnumLoops)	\
    (This)->lpVtbl -> Clone(This,pEnumLoops)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMLoops_RemoteNext_Proxy( 
    IEnumDMLoops __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [length_is][size_is][out] */ LPDMLOOP __RPC_FAR *pLoop,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMLoops_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMLoops_Skip_Proxy( 
    IEnumDMLoops __RPC_FAR * This,
    /* [in] */ ULONG cLoops);


void __RPC_STUB IEnumDMLoops_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMLoops_Reset_Proxy( 
    IEnumDMLoops __RPC_FAR * This);


void __RPC_STUB IEnumDMLoops_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMLoops_Clone_Proxy( 
    IEnumDMLoops __RPC_FAR * This,
    /* [out] */ LPENUM_DMLOOPS __RPC_FAR *pEnumLoops);


void __RPC_STUB IEnumDMLoops_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMLoops_INTERFACE_DEFINED__ */


#ifndef __IEnumDMEdgeUses_INTERFACE_DEFINED__
#define __IEnumDMEdgeUses_INTERFACE_DEFINED__

/* interface IEnumDMEdgeUses */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMEdgeUses;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D227-0000-0000-C000-000000000046")
    IEnumDMEdgeUses : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cEdgeUse,
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cEdgeUses) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMEdgeUsesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMEdgeUses __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMEdgeUses __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMEdgeUses __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMEdgeUses __RPC_FAR * This,
            /* [in] */ ULONG cEdgeUse,
            /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMEdgeUses __RPC_FAR * This,
            /* [in] */ ULONG cEdgeUses);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMEdgeUses __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMEdgeUses __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses);
        
        END_INTERFACE
    } IEnumDMEdgeUsesVtbl;

    interface IEnumDMEdgeUses
    {
        CONST_VTBL struct IEnumDMEdgeUsesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMEdgeUses_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMEdgeUses_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMEdgeUses_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMEdgeUses_Next(This,cEdgeUse,pEdgeUse,pcFetched)	\
    (This)->lpVtbl -> Next(This,cEdgeUse,pEdgeUse,pcFetched)

#define IEnumDMEdgeUses_Skip(This,cEdgeUses)	\
    (This)->lpVtbl -> Skip(This,cEdgeUses)

#define IEnumDMEdgeUses_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMEdgeUses_Clone(This,pEnumEdgeUses)	\
    (This)->lpVtbl -> Clone(This,pEnumEdgeUses)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_RemoteNext_Proxy( 
    IEnumDMEdgeUses __RPC_FAR * This,
    /* [in] */ ULONG cEdgeUse,
    /* [length_is][size_is][out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMEdgeUses_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_Skip_Proxy( 
    IEnumDMEdgeUses __RPC_FAR * This,
    /* [in] */ ULONG cEdgeUses);


void __RPC_STUB IEnumDMEdgeUses_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_Reset_Proxy( 
    IEnumDMEdgeUses __RPC_FAR * This);


void __RPC_STUB IEnumDMEdgeUses_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_Clone_Proxy( 
    IEnumDMEdgeUses __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGEUSES __RPC_FAR *pEnumEdgeUses);


void __RPC_STUB IEnumDMEdgeUses_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMEdgeUses_INTERFACE_DEFINED__ */


#ifndef __IEnumDMEdges_INTERFACE_DEFINED__
#define __IEnumDMEdges_INTERFACE_DEFINED__

/* interface IEnumDMEdges */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMEdges;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D228-0000-0000-C000-000000000046")
    IEnumDMEdges : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cEdge,
            /* [out] */ LPDMEDGE __RPC_FAR *pEdge,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cEdges) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMEdgesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMEdges __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMEdges __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMEdges __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMEdges __RPC_FAR * This,
            /* [in] */ ULONG cEdge,
            /* [out] */ LPDMEDGE __RPC_FAR *pEdge,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMEdges __RPC_FAR * This,
            /* [in] */ ULONG cEdges);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMEdges __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMEdges __RPC_FAR * This,
            /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);
        
        END_INTERFACE
    } IEnumDMEdgesVtbl;

    interface IEnumDMEdges
    {
        CONST_VTBL struct IEnumDMEdgesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMEdges_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMEdges_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMEdges_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMEdges_Next(This,cEdge,pEdge,pcFetched)	\
    (This)->lpVtbl -> Next(This,cEdge,pEdge,pcFetched)

#define IEnumDMEdges_Skip(This,cEdges)	\
    (This)->lpVtbl -> Skip(This,cEdges)

#define IEnumDMEdges_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMEdges_Clone(This,pEnumEdges)	\
    (This)->lpVtbl -> Clone(This,pEnumEdges)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMEdges_RemoteNext_Proxy( 
    IEnumDMEdges __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [length_is][size_is][out] */ LPDMEDGE __RPC_FAR *pEdge,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMEdges_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdges_Skip_Proxy( 
    IEnumDMEdges __RPC_FAR * This,
    /* [in] */ ULONG cEdges);


void __RPC_STUB IEnumDMEdges_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdges_Reset_Proxy( 
    IEnumDMEdges __RPC_FAR * This);


void __RPC_STUB IEnumDMEdges_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMEdges_Clone_Proxy( 
    IEnumDMEdges __RPC_FAR * This,
    /* [out] */ LPENUM_DMEDGES __RPC_FAR *pEnumEdges);


void __RPC_STUB IEnumDMEdges_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMEdges_INTERFACE_DEFINED__ */


#ifndef __IEnumDMVertices_INTERFACE_DEFINED__
#define __IEnumDMVertices_INTERFACE_DEFINED__

/* interface IEnumDMVertices */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMVertices;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D229-0000-0000-C000-000000000046")
    IEnumDMVertices : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cVertices,
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cVertices) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMVerticesVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMVertices __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMVertices __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMVertices __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMVertices __RPC_FAR * This,
            /* [in] */ ULONG cVertices,
            /* [out] */ LPDMVERTEX __RPC_FAR *pVertex,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMVertices __RPC_FAR * This,
            /* [in] */ ULONG cVertices);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMVertices __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMVertices __RPC_FAR * This,
            /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices);
        
        END_INTERFACE
    } IEnumDMVerticesVtbl;

    interface IEnumDMVertices
    {
        CONST_VTBL struct IEnumDMVerticesVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMVertices_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMVertices_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMVertices_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMVertices_Next(This,cVertices,pVertex,pcFetched)	\
    (This)->lpVtbl -> Next(This,cVertices,pVertex,pcFetched)

#define IEnumDMVertices_Skip(This,cVertices)	\
    (This)->lpVtbl -> Skip(This,cVertices)

#define IEnumDMVertices_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMVertices_Clone(This,pEnumVertices)	\
    (This)->lpVtbl -> Clone(This,pEnumVertices)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMVertices_RemoteNext_Proxy( 
    IEnumDMVertices __RPC_FAR * This,
    /* [in] */ ULONG cVertices,
    /* [length_is][size_is][out] */ LPDMVERTEX __RPC_FAR *pVertex,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMVertices_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMVertices_Skip_Proxy( 
    IEnumDMVertices __RPC_FAR * This,
    /* [in] */ ULONG cVertices);


void __RPC_STUB IEnumDMVertices_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMVertices_Reset_Proxy( 
    IEnumDMVertices __RPC_FAR * This);


void __RPC_STUB IEnumDMVertices_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMVertices_Clone_Proxy( 
    IEnumDMVertices __RPC_FAR * This,
    /* [out] */ LPENUM_DMVERTICES __RPC_FAR *pEnumVertices);


void __RPC_STUB IEnumDMVertices_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMVertices_INTERFACE_DEFINED__ */


#ifndef __IDMAlternateSurfaceBody_INTERFACE_DEFINED__
#define __IDMAlternateSurfaceBody_INTERFACE_DEFINED__

/* interface IDMAlternateSurfaceBody */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMAlternateSurfaceBody;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D22A-0000-0000-C000-000000000046")
    IDMAlternateSurfaceBody : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetAlternateBody( 
            /* [in] */ DWORD alternateForm,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *ppAltBody) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMAlternateSurfaceBodyVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMAlternateSurfaceBody __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMAlternateSurfaceBody __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMAlternateSurfaceBody __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetAlternateBody )( 
            IDMAlternateSurfaceBody __RPC_FAR * This,
            /* [in] */ DWORD alternateForm,
            /* [out] */ LPDMSURFACEBODY __RPC_FAR *ppAltBody);
        
        END_INTERFACE
    } IDMAlternateSurfaceBodyVtbl;

    interface IDMAlternateSurfaceBody
    {
        CONST_VTBL struct IDMAlternateSurfaceBodyVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMAlternateSurfaceBody_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMAlternateSurfaceBody_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMAlternateSurfaceBody_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMAlternateSurfaceBody_GetAlternateBody(This,alternateForm,ppAltBody)	\
    (This)->lpVtbl -> GetAlternateBody(This,alternateForm,ppAltBody)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMAlternateSurfaceBody_GetAlternateBody_Proxy( 
    IDMAlternateSurfaceBody __RPC_FAR * This,
    /* [in] */ DWORD alternateForm,
    /* [out] */ LPDMSURFACEBODY __RPC_FAR *ppAltBody);


void __RPC_STUB IDMAlternateSurfaceBody_GetAlternateBody_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMAlternateSurfaceBody_INTERFACE_DEFINED__ */


/* interface __MIDL_itf_gtfordm_0109 */
/* [local] */ 

DEFINE_GUID(IID_IDMReferenceKey, 0x0002D22B, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IDMReference, 0x0002D22C, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
DEFINE_GUID(IID_IEnumDMReferenceKeys, 0x0002D22D, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
typedef /* [unique] */ IDMReferenceKey __RPC_FAR *LPDMREFKEY;

typedef /* [unique] */ IDMReference __RPC_FAR *LPDMREF;

typedef /* [unique] */ IEnumDMReferenceKeys __RPC_FAR *LPENUM_DMREFKEYS;



extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0109_v0_0_c_ifspec;
extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0109_v0_0_s_ifspec;

#ifndef __IDMReferenceKey_INTERFACE_DEFINED__
#define __IDMReferenceKey_INTERFACE_DEFINED__

/* interface IDMReferenceKey */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMReferenceKey;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D22B-0000-0000-C000-000000000046")
    IDMReferenceKey : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE GetObjectType( 
            /* [out] */ IID __RPC_FAR *pRefIID) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetKeySize( 
            /* [out] */ ULONG __RPC_FAR *pcbKeySize) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetKey( 
            /* [in] */ ULONG cbKeySize,
            /* [out][size_is] */ BYTE __RPC_FAR *pKey) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetModifyCounter( 
            /* [out] */ long __RPC_FAR *pCounter) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMReferenceKeyVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMReferenceKey __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMReferenceKey __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMReferenceKey __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetObjectType )( 
            IDMReferenceKey __RPC_FAR * This,
            /* [out] */ IID __RPC_FAR *pRefIID);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetKeySize )( 
            IDMReferenceKey __RPC_FAR * This,
            /* [out] */ ULONG __RPC_FAR *pcbKeySize);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetKey )( 
            IDMReferenceKey __RPC_FAR * This,
            /* [in] */ ULONG cbKeySize,
            /* [out][size_is] */ BYTE __RPC_FAR *pKey);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetModifyCounter )( 
            IDMReferenceKey __RPC_FAR * This,
            /* [out] */ long __RPC_FAR *pCounter);
        
        END_INTERFACE
    } IDMReferenceKeyVtbl;

    interface IDMReferenceKey
    {
        CONST_VTBL struct IDMReferenceKeyVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMReferenceKey_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMReferenceKey_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMReferenceKey_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMReferenceKey_GetObjectType(This,pRefIID)	\
    (This)->lpVtbl -> GetObjectType(This,pRefIID)

#define IDMReferenceKey_GetKeySize(This,pcbKeySize)	\
    (This)->lpVtbl -> GetKeySize(This,pcbKeySize)

#define IDMReferenceKey_GetKey(This,cbKeySize,pKey)	\
    (This)->lpVtbl -> GetKey(This,cbKeySize,pKey)

#define IDMReferenceKey_GetModifyCounter(This,pCounter)	\
    (This)->lpVtbl -> GetModifyCounter(This,pCounter)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMReferenceKey_GetObjectType_Proxy( 
    IDMReferenceKey __RPC_FAR * This,
    /* [out] */ IID __RPC_FAR *pRefIID);


void __RPC_STUB IDMReferenceKey_GetObjectType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMReferenceKey_GetKeySize_Proxy( 
    IDMReferenceKey __RPC_FAR * This,
    /* [out] */ ULONG __RPC_FAR *pcbKeySize);


void __RPC_STUB IDMReferenceKey_GetKeySize_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMReferenceKey_GetKey_Proxy( 
    IDMReferenceKey __RPC_FAR * This,
    /* [in] */ ULONG cbKeySize,
    /* [out][size_is] */ BYTE __RPC_FAR *pKey);


void __RPC_STUB IDMReferenceKey_GetKey_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMReferenceKey_GetModifyCounter_Proxy( 
    IDMReferenceKey __RPC_FAR * This,
    /* [out] */ long __RPC_FAR *pCounter);


void __RPC_STUB IDMReferenceKey_GetModifyCounter_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMReferenceKey_INTERFACE_DEFINED__ */


#ifndef __IDMReference_INTERFACE_DEFINED__
#define __IDMReference_INTERFACE_DEFINED__

/* interface IDMReference */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMReference;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D22C-0000-0000-C000-000000000046")
    IDMReference : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE BindKeyToInterface( 
            /* [in] */ REFIID IDMInterface,
            /* [in] */ DWORD cSizeInBytes,
            /* [size_is][in] */ BYTE __RPC_FAR *pKey,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMReferenceVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMReference __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMReference __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMReference __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *BindKeyToInterface )( 
            IDMReference __RPC_FAR * This,
            /* [in] */ REFIID IDMInterface,
            /* [in] */ DWORD cSizeInBytes,
            /* [size_is][in] */ BYTE __RPC_FAR *pKey,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        END_INTERFACE
    } IDMReferenceVtbl;

    interface IDMReference
    {
        CONST_VTBL struct IDMReferenceVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMReference_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMReference_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMReference_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMReference_BindKeyToInterface(This,IDMInterface,cSizeInBytes,pKey,ppvObject)	\
    (This)->lpVtbl -> BindKeyToInterface(This,IDMInterface,cSizeInBytes,pKey,ppvObject)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMReference_BindKeyToInterface_Proxy( 
    IDMReference __RPC_FAR * This,
    /* [in] */ REFIID IDMInterface,
    /* [in] */ DWORD cSizeInBytes,
    /* [size_is][in] */ BYTE __RPC_FAR *pKey,
    /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);


void __RPC_STUB IDMReference_BindKeyToInterface_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMReference_INTERFACE_DEFINED__ */


#ifndef __IEnumDMReferenceKeys_INTERFACE_DEFINED__
#define __IEnumDMReferenceKeys_INTERFACE_DEFINED__

/* interface IEnumDMReferenceKeys */
/* [object][uuid] */ 


EXTERN_C const IID IID_IEnumDMReferenceKeys;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D22D-0000-0000-C000-000000000046")
    IEnumDMReferenceKeys : public IUnknown
    {
    public:
        virtual /* [local] */ HRESULT STDMETHODCALLTYPE Next( 
            /* [in] */ ULONG cRefKey,
            /* [out] */ LPDMREFKEY __RPC_FAR *pRefKey,
            /* [out] */ ULONG __RPC_FAR *pcFetched) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Skip( 
            /* [in] */ ULONG cRefKey) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Reset( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE Clone( 
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKeys) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IEnumDMReferenceKeysVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IEnumDMReferenceKeys __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IEnumDMReferenceKeys __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IEnumDMReferenceKeys __RPC_FAR * This);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Next )( 
            IEnumDMReferenceKeys __RPC_FAR * This,
            /* [in] */ ULONG cRefKey,
            /* [out] */ LPDMREFKEY __RPC_FAR *pRefKey,
            /* [out] */ ULONG __RPC_FAR *pcFetched);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Skip )( 
            IEnumDMReferenceKeys __RPC_FAR * This,
            /* [in] */ ULONG cRefKey);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Reset )( 
            IEnumDMReferenceKeys __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Clone )( 
            IEnumDMReferenceKeys __RPC_FAR * This,
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKeys);
        
        END_INTERFACE
    } IEnumDMReferenceKeysVtbl;

    interface IEnumDMReferenceKeys
    {
        CONST_VTBL struct IEnumDMReferenceKeysVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEnumDMReferenceKeys_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IEnumDMReferenceKeys_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IEnumDMReferenceKeys_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IEnumDMReferenceKeys_Next(This,cRefKey,pRefKey,pcFetched)	\
    (This)->lpVtbl -> Next(This,cRefKey,pRefKey,pcFetched)

#define IEnumDMReferenceKeys_Skip(This,cRefKey)	\
    (This)->lpVtbl -> Skip(This,cRefKey)

#define IEnumDMReferenceKeys_Reset(This)	\
    (This)->lpVtbl -> Reset(This)

#define IEnumDMReferenceKeys_Clone(This,pEnumReferenceKeys)	\
    (This)->lpVtbl -> Clone(This,pEnumReferenceKeys)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_RemoteNext_Proxy( 
    IEnumDMReferenceKeys __RPC_FAR * This,
    /* [in] */ ULONG cRefKey,
    /* [length_is][size_is][out] */ LPDMREFKEY __RPC_FAR *pRefKey,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


void __RPC_STUB IEnumDMReferenceKeys_RemoteNext_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_Skip_Proxy( 
    IEnumDMReferenceKeys __RPC_FAR * This,
    /* [in] */ ULONG cRefKey);


void __RPC_STUB IEnumDMReferenceKeys_Skip_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_Reset_Proxy( 
    IEnumDMReferenceKeys __RPC_FAR * This);


void __RPC_STUB IEnumDMReferenceKeys_Reset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_Clone_Proxy( 
    IEnumDMReferenceKeys __RPC_FAR * This,
    /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKeys);


void __RPC_STUB IEnumDMReferenceKeys_Clone_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IEnumDMReferenceKeys_INTERFACE_DEFINED__ */


/* interface __MIDL_itf_gtfordm_0112 */
/* [local] */ 

DEFINE_GUID(IID_IDMGeometricLocate, 0x0002D22E, 0x0000, 0x0000, 0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
typedef /* [unique] */ IDMGeometricLocate __RPC_FAR *LPDMGEOMETRICLOCATE;

typedef struct  tagDMBoreLine
    {
    double m_point[ 3 ];
    double m_direction[ 3 ];
    double m_front;
    double m_back;
    double m_radius;
    }	DMBORELINE;

typedef DMBORELINE __RPC_FAR *LPDMBORELINE;

typedef 
enum tagDMSELECTTYPE
    {	SELECTTYPE_INSIDE	= 0x1,
	SELECTTYPE_OVERLAP	= 0x2
    }	DMSELECTTYPE;

typedef struct  tagDMShape
    {
    ULONG m_nPoints;
    double __RPC_FAR *m_pPoints;
    double m_direction[ 3 ];
    double m_front;
    double m_back;
    DWORD m_type;
    }	DMSHAPE;

typedef DMSHAPE __RPC_FAR *LPDMSHAPE;



extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0112_v0_0_c_ifspec;
extern RPC_IF_HANDLE __MIDL_itf_gtfordm_0112_v0_0_s_ifspec;

#ifndef __IDMGeometricLocate_INTERFACE_DEFINED__
#define __IDMGeometricLocate_INTERFACE_DEFINED__

/* interface IDMGeometricLocate */
/* [object][uuid] */ 


EXTERN_C const IID IID_IDMGeometricLocate;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("0002D22E-0000-0000-C000-000000000046")
    IDMGeometricLocate : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE PointLocate( 
            /* [in] */ LPDMBORELINE pBoreline,
            /* [in] */ ULONG nTypes,
            /* [in][size_is] */ IID __RPC_FAR types[  ],
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE ShapeLocate( 
            /* [in] */ LPDMSHAPE pShape,
            /* [in] */ ULONG nTypes,
            /* [in][size_is] */ IID __RPC_FAR types[  ],
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDMGeometricLocateVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDMGeometricLocate __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDMGeometricLocate __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDMGeometricLocate __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PointLocate )( 
            IDMGeometricLocate __RPC_FAR * This,
            /* [in] */ LPDMBORELINE pBoreline,
            /* [in] */ ULONG nTypes,
            /* [in][size_is] */ IID __RPC_FAR types[  ],
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ShapeLocate )( 
            IDMGeometricLocate __RPC_FAR * This,
            /* [in] */ LPDMSHAPE pShape,
            /* [in] */ ULONG nTypes,
            /* [in][size_is] */ IID __RPC_FAR types[  ],
            /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey);
        
        END_INTERFACE
    } IDMGeometricLocateVtbl;

    interface IDMGeometricLocate
    {
        CONST_VTBL struct IDMGeometricLocateVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDMGeometricLocate_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDMGeometricLocate_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDMGeometricLocate_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDMGeometricLocate_PointLocate(This,pBoreline,nTypes,types,pEnumReferenceKey)	\
    (This)->lpVtbl -> PointLocate(This,pBoreline,nTypes,types,pEnumReferenceKey)

#define IDMGeometricLocate_ShapeLocate(This,pShape,nTypes,types,pEnumReferenceKey)	\
    (This)->lpVtbl -> ShapeLocate(This,pShape,nTypes,types,pEnumReferenceKey)

#endif /* COBJMACROS */


#endif 	/* C style interface */



HRESULT STDMETHODCALLTYPE IDMGeometricLocate_PointLocate_Proxy( 
    IDMGeometricLocate __RPC_FAR * This,
    /* [in] */ LPDMBORELINE pBoreline,
    /* [in] */ ULONG nTypes,
    /* [in][size_is] */ IID __RPC_FAR types[  ],
    /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey);


void __RPC_STUB IDMGeometricLocate_PointLocate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


HRESULT STDMETHODCALLTYPE IDMGeometricLocate_ShapeLocate_Proxy( 
    IDMGeometricLocate __RPC_FAR * This,
    /* [in] */ LPDMSHAPE pShape,
    /* [in] */ ULONG nTypes,
    /* [in][size_is] */ IID __RPC_FAR types[  ],
    /* [out] */ LPENUM_DMREFKEYS __RPC_FAR *pEnumReferenceKey);


void __RPC_STUB IDMGeometricLocate_ShapeLocate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDMGeometricLocate_INTERFACE_DEFINED__ */


/* Additional Prototypes for ALL interfaces */

/* [local] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetParamAtPoint_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [in] */ double __RPC_FAR *pPoints,
    /* [in] */ double __RPC_FAR *pGuessParams,
    /* [out] */ double __RPC_FAR *pMaxDeviations,
    /* [out] */ double __RPC_FAR *pParams,
    /* [out] */ DWORD __RPC_FAR *pFlags);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetParamAtPoint_Stub( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetTangents_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pUTangents,
    /* [out] */ double __RPC_FAR *pVTangents);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetTangents_Stub( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nUTangents,
    /* [size_is][out] */ double __RPC_FAR *pUTangents,
    /* [in] */ ULONG nVTangents,
    /* [size_is][out] */ double __RPC_FAR *pVTangents);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetCurvatures_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pMaxTangents,
    /* [out] */ double __RPC_FAR *pMaxCurvatures,
    /* [out] */ double __RPC_FAR *pMinCurvatures);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetCurvatures_Stub( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nMaxTangents,
    /* [size_is][out] */ double __RPC_FAR *pMaxTangents,
    /* [in] */ ULONG nMaxCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pMaxCurvatures,
    /* [in] */ ULONG nMinCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pMinCurvatures);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetDerivatives_Proxy( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pUPartials,
    /* [out] */ double __RPC_FAR *pVPartials,
    /* [out] */ double __RPC_FAR *pUUPartials,
    /* [out] */ double __RPC_FAR *pUVPartials,
    /* [out] */ double __RPC_FAR *pVVPartials,
    /* [out] */ double __RPC_FAR *pUUUPartials,
    /* [out] */ double __RPC_FAR *pVVVPartials);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMSurface_GetDerivatives_Stub( 
    IDMSurface __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUPartials,
    /* [in] */ ULONG nVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVPartials,
    /* [in] */ ULONG nUUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUUPartials,
    /* [in] */ ULONG nUVPartials,
    /* [size_is][out] */ double __RPC_FAR *pUVPartials,
    /* [in] */ ULONG nVVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVVPartials,
    /* [in] */ ULONG nUUUPartials,
    /* [size_is][out] */ double __RPC_FAR *pUUUPartials,
    /* [in] */ ULONG nVVVPartials,
    /* [size_is][out] */ double __RPC_FAR *pVVVPartials);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetParamAtPoint_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [in] */ double __RPC_FAR *pPoints,
    /* [in] */ double __RPC_FAR *pGuessParams,
    /* [out] */ double __RPC_FAR *pMaxDeviations,
    /* [out] */ double __RPC_FAR *pParams,
    /* [out] */ DWORD __RPC_FAR *pFlags);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetParamAtPoint_Stub( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetCurvature_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pDirections,
    /* [out] */ double __RPC_FAR *pCurvatures);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetCurvature_Stub( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nDirections,
    /* [size_is][out] */ double __RPC_FAR *pDirections,
    /* [in] */ ULONG nCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pCurvatures);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetDerivatives_Proxy( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pFirstDerivs,
    /* [out] */ double __RPC_FAR *pSecondDerivs,
    /* [out] */ double __RPC_FAR *pThirdDerivs);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve_GetDerivatives_Stub( 
    IDMCurve __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFirstDerivs,
    /* [size_is][out] */ double __RPC_FAR *pFirstDerivs,
    /* [in] */ ULONG nSecondDerivs,
    /* [size_is][out] */ double __RPC_FAR *pSecondDerivs,
    /* [in] */ ULONG nThirdDerivs,
    /* [size_is][out] */ double __RPC_FAR *pThirdDerivs);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetParamAtPoint_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [in] */ double __RPC_FAR *pPoints,
    /* [in] */ double __RPC_FAR *pGuessParams,
    /* [out] */ double __RPC_FAR *pMaxDeviations,
    /* [out] */ double __RPC_FAR *pParams,
    /* [out] */ DWORD __RPC_FAR *pFlags);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetParamAtPoint_Stub( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nPoints,
    /* [size_is][in] */ double __RPC_FAR *pPoints,
    /* [in] */ ULONG nGuessParams,
    /* [size_is][in] */ double __RPC_FAR *pGuessParams,
    /* [in] */ ULONG nMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pMaxDeviations,
    /* [size_is][out] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFlags,
    /* [size_is][out] */ DWORD __RPC_FAR *pFlags);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetCurvature_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pDirections,
    /* [out] */ double __RPC_FAR *pCurvatures);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetCurvature_Stub( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nDirections,
    /* [size_is][out] */ double __RPC_FAR *pDirections,
    /* [in] */ ULONG nCurvatures,
    /* [size_is][out] */ double __RPC_FAR *pCurvatures);

/* [local] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetDerivatives_Proxy( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [in] */ double __RPC_FAR *pParams,
    /* [out] */ double __RPC_FAR *pFirstDerivs,
    /* [out] */ double __RPC_FAR *pSecondDerivs,
    /* [out] */ double __RPC_FAR *pThirdDerivs);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IDMCurve2D_GetDerivatives_Stub( 
    IDMCurve2D __RPC_FAR * This,
    /* [in] */ ULONG nParams,
    /* [size_is][in] */ double __RPC_FAR *pParams,
    /* [in] */ ULONG nFirstDerivs,
    /* [size_is][out] */ double __RPC_FAR *pFirstDerivs,
    /* [in] */ ULONG nSecondDerivs,
    /* [size_is][out] */ double __RPC_FAR *pSecondDerivs,
    /* [in] */ ULONG nThirdDerivs,
    /* [size_is][out] */ double __RPC_FAR *pThirdDerivs);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_Next_Proxy( 
    IEnumDMSurfaceBodies __RPC_FAR * This,
    /* [in] */ ULONG cSurfaceBody,
    /* [out] */ LPDMSURFACEBODY __RPC_FAR *pSurfaceBody,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMSurfaceBodies_Next_Stub( 
    IEnumDMSurfaceBodies __RPC_FAR * This,
    /* [in] */ ULONG cSurfaceBody,
    /* [length_is][size_is][out] */ LPDMSURFACEBODY __RPC_FAR *pSurfaceBody,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMShells_Next_Proxy( 
    IEnumDMShells __RPC_FAR * This,
    /* [in] */ ULONG cShells,
    /* [out] */ LPDMSHELL __RPC_FAR *pShell,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMShells_Next_Stub( 
    IEnumDMShells __RPC_FAR * This,
    /* [in] */ ULONG cShells,
    /* [length_is][size_is][out] */ LPDMSHELL __RPC_FAR *pShell,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMFaces_Next_Proxy( 
    IEnumDMFaces __RPC_FAR * This,
    /* [in] */ ULONG cFaces,
    /* [out] */ LPDMFACE __RPC_FAR *pFace,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMFaces_Next_Stub( 
    IEnumDMFaces __RPC_FAR * This,
    /* [in] */ ULONG cFaces,
    /* [length_is][size_is][out] */ LPDMFACE __RPC_FAR *pFace,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMLoops_Next_Proxy( 
    IEnumDMLoops __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [out] */ LPDMLOOP __RPC_FAR *pLoop,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMLoops_Next_Stub( 
    IEnumDMLoops __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [length_is][size_is][out] */ LPDMLOOP __RPC_FAR *pLoop,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_Next_Proxy( 
    IEnumDMEdgeUses __RPC_FAR * This,
    /* [in] */ ULONG cEdgeUse,
    /* [out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMEdgeUses_Next_Stub( 
    IEnumDMEdgeUses __RPC_FAR * This,
    /* [in] */ ULONG cEdgeUse,
    /* [length_is][size_is][out] */ LPDMEDGEUSE __RPC_FAR *pEdgeUse,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMEdges_Next_Proxy( 
    IEnumDMEdges __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [out] */ LPDMEDGE __RPC_FAR *pEdge,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMEdges_Next_Stub( 
    IEnumDMEdges __RPC_FAR * This,
    /* [in] */ ULONG cEdge,
    /* [length_is][size_is][out] */ LPDMEDGE __RPC_FAR *pEdge,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMVertices_Next_Proxy( 
    IEnumDMVertices __RPC_FAR * This,
    /* [in] */ ULONG cVertices,
    /* [out] */ LPDMVERTEX __RPC_FAR *pVertex,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMVertices_Next_Stub( 
    IEnumDMVertices __RPC_FAR * This,
    /* [in] */ ULONG cVertices,
    /* [length_is][size_is][out] */ LPDMVERTEX __RPC_FAR *pVertex,
    /* [out] */ ULONG __RPC_FAR *pcFetched);

/* [local] */ HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_Next_Proxy( 
    IEnumDMReferenceKeys __RPC_FAR * This,
    /* [in] */ ULONG cRefKey,
    /* [out] */ LPDMREFKEY __RPC_FAR *pRefKey,
    /* [out] */ ULONG __RPC_FAR *pcFetched);


/* [call_as] */ HRESULT STDMETHODCALLTYPE IEnumDMReferenceKeys_Next_Stub( 
    IEnumDMReferenceKeys __RPC_FAR * This,
    /* [in] */ ULONG cRefKey,
    /* [length_is][size_is][out] */ LPDMREFKEY __RPC_FAR *pRefKey,
    /* [out] */ ULONG __RPC_FAR *pcFetched);



/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif
