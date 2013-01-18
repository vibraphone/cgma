//-------------------------------------------------------------------------
// Filename      : SplitSurfaceVirtual.hpp
//
// Purpose       : Split surface(s) and create virtual surfaces
//
// Special Notes : 
//
// Creator       : KGM
//
// Creation Date : 30-Jun-2005
//-------------------------------------------------------------------------

#ifndef SPLITSURFACEVIRTUAL_HPP
#define SPLITSURFACEVIRTUAL_HPP

class RefVertex;
class RefEdge;
class RefFace;
class CubitVector;
template <class X> class DLIList;

#include "CubitDefines.h"

namespace SplitSurfaceVirtual
{

  CubitStatus split_surface_virtual( RefFace *ref_face_ptr,
                                     DLIList<CubitVector*> &locations,
                                     DLIList<DLIList<CubitVector*>*> &vec_lists );
    //- Split a surface and create virtual surfaces as a result.

  CubitStatus split_surfaces_virtual( DLIList<RefFace*> &ref_face_list, 
                                       int num_segs, double fraction,
                                       double distance, RefEdge *from_curve_ptr,
                                       DLIList<RefVertex*> &corner_vertex_list,
                                       DLIList<RefVertex*> &through_vertex_list,
                                       RefEdge *curve_dir_ptr, CubitBoolean preview_flg,
                                       CubitBoolean create_ref_edges_flg );
};

#endif 

