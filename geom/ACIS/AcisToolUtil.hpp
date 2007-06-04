//-------------------------------------------------------------------------
// Filename      : AcisToolUtil.hpp
//
// Purpose       : Provide old RefEntity-based functions that were removed
//                 from Acis{Query|Modify}Engine for use by misc. Acis*Tool
//                 classes that still have a RefEntity-based interface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/23/03
//-------------------------------------------------------------------------


#ifndef ACIS_TOOL_UTIL_HPP
#define ACIS_TOOL_UTIL_HPP

#include "CubitDefines.h"  /* for CubitStatus */

template <class X> class DLIList;

class Body;
class RefEdge;
class RefFace;

class ENTITY;
class BODY;
class FACE;
class EDGE;

class CurveACIS;
class SurfaceACIS;

class AcisToolUtil
{
  public:
  
    static Body* get_new_Body( Body* old_Body, 
                               BODY* old_BODY, 
                               BODY* new_body,
                               bool keep_old,
                               bool topo_check = false,
                               bool delete_old = true );
    
    static CubitStatus 
    get_ACIS_Curves( DLIList<RefEdge*>& ref_edge_list,
                     DLIList<CurveACIS*>& curve_list );
                                 
    static CubitStatus
    get_ACIS_Surfaces( DLIList<RefFace*>& face_list,
                       DLIList<SurfaceACIS*>& surface_list );
    
#ifdef BOYD14
    static CubitStatus 
    get_copied_EDGES_of_body( DLIList<RefEdge*>& ref_edge_list,
                              DLIList<EDGE*>& copied_EDGEs,
                              BODY*& copied_body_ptr );
#endif
                                          
    static CubitStatus 
    get_copied_FACES_of_body( DLIList<RefFace*>& ref_face_list,
                              DLIList<FACE*>& copied_FACEs,
                              BODY*& copied_body_ptr );
    
    static CubitStatus 
    get_copied_FACES_of_body( DLIList<RefFace*>& ref_face_list,
                              DLIList<FACE*>& copied_FACEs,
                              DLIList<RefFace*>& removed_face_list,
                              BODY*& copied_body_ptr );
    
    static Body* get_body_of_ENTITY(ENTITY*);

};

#endif
