
#ifndef _REMOVE_BLENDS_
#define _REMOVE_BLENDS_

#include "CubitUtil.hpp"
#include "DLIList.hpp"

class RefVertex;
class RefEdge;
class RefFace;
class CubitVector;

namespace RemoveBlends
{

    CubitStatus remove_blend(RefFace* ref_face,
                             DLIList<CubitVector*>& locations,
                             DLIList<RefFace*>& composite_faces);

    CubitStatus remove_blends(DLIList<RefFace*>& ref_face_list,
                              int num_segs, double fraction, double distance,
                              RefEdge* from_curve_ptr,
                              DLIList<RefVertex*>& corner_vertex_list,
                              DLIList<RefVertex*>& through_vertex_list,
                              RefEdge *curve_dir_ptr,
                              CubitBoolean preview_flg,
                              DLIList<CubitVector*>& locations,
                              DLIList<RefFace*>& composite_faces);
};
#endif

