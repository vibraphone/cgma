//-Class:       BoundaryContstrainTool
//-Description: recover edges from a triangle mesh.  Currently 2D meshes only
//-Owner:       Steve Owen
//-Altered by: Jit Ken Tan (07/08/02)
//-Checked by:

#ifndef BOUNDARY_CONSTRAIN_TOOL_HPP
#define BOUNDARY_CONSTRAIN_TOOL_HPP

#include "DLIList.hpp"
class CubitVector;

typedef enum SwapStatus { SWAP_FAILURE, SWAP_SUCCESS, SWAP_INVALID } SwapStatus;
typedef enum IntersectionType { AT_BEGIN, AT_EDGE, AT_MID, AT_END, NO_ISECT } IntersectionType;

template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> class BoundaryConstrainTool
{
private:
  DLIList <TRI *> *facetList;
  SURF *refFacePtr;
  double zeroTol;
  DLIList<EDGE *> edgeCrossQueue;

  CubitStatus get_crossing_edges( CubitVector &edge_vec, 
                                  NODE *n0_ptr, 
                                  NODE *n1_ptr);
   //- generate the initial list of all edges that cross the edge_vec
   //- between n0 and n1

  IntersectionType intersect_from_node(CubitVector &edge_vec, 
                                       NODE *n0_ptr, 
                                       TRI *&tri_ptr, 
                                       EDGE *&edge_ptr, 
                                       NODE *&node_ptr );
   //- determine next intersection assuming start from n0

  IntersectionType intersect_from_edge( CubitVector &edge_vec, 
                                        NODE *n0_ptr, 
                                        NODE *n1_ptr,
                                        TRI *&tri_ptr, 
                                        EDGE *&edge_ptr, 
                                        NODE *&node_ptr );
   //- determine next intersection assuming start from edge_ptr

  CubitStatus node_at_mid( CubitVector &edge_vec, 
                           TRI *tri_ptr, 
                           NODE *node_ptr );
   //- he vector hit a node on its way to the n1

  SwapStatus swap_edge( EDGE *&edge_ptr );
   //- attempt to swap a single edge.  Check for valid swap 
   //- before doing so

  CubitStatus edge_intersected( CubitVector &edge_vec, 
                                EDGE *edge_ptr, 
                                NODE *n0_ptr, 
                                NODE *n1_ptr );
   //- check to see if the edge we just generated crosses the
   //- edge_vec.  If so, then add it to the queue 

public:
  BoundaryConstrainTool();
  BoundaryConstrainTool(SURF *ref_face_ptr);
  ~BoundaryConstrainTool();
    //- constructor and desructor

  CubitStatus recover_edge(NODE *n0_ptr,NODE *n1_ptr,
                           EDGE *&recovered_edge_ptr,
			   DLIList <TRI *> *facet_list);

    //- flip edges in the triangulation to recover the edge between
    //- n0_ptr and n1_ptr. Assumes 2D (x-y-0)
};

#include "BoundaryConstrainTool.cpp"

#endif

