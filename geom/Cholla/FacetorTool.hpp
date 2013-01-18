//- Class: FacetorTool
//- Description: FacetorTool is a template class that accepts both CubitNode, CubitTri, CubitEdge classes 
//-              and CubitPoint, CubitFacet, CubitFacetEdge classes.  CubitPoint, CubitFacet, and 
//-              CubitFacetEdge classes have been modified to account for the difference between the mdb
//=              entities and geometric entities.
//-------------------------------------------------------------------------
// Filename      : FacetorTool.hpp
//
// Purpose       : 2D Delaunay Mesher
//
// Creator       : Christopher Hynes
//
// Creation Date : 5/31/2002
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

#ifndef FACETOR_TOOL_HPP
#define FACETOR_TOOL_HPP

#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"

class CubitFacet;
class CubitFacetEdge;
class CubitPoint;
class ParamTool;

//template <class X> class DLIList;
#include "DLIList.hpp"

template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> 
class FacetorTool
{
public:
  FacetorTool(){};

  FacetorTool( SURF *ref_face_ptr,
	       DLIList<NODE *> &bounding_nodes,
	       NODE **boundary_edge_start_nodes, 
	       NODE **boundary_edge_end_nodes, 
	       int num_boundary_edges,
	       SIZEFUNC *sizing_function,
		   ParamTool *ptool);

  ~FacetorTool(void);

   CubitStatus mesh_surf(DLIList<TRI *> &facet_list);
   //- mesh the surface
   //- returns a facet_list of the mesh entities.  If using mesh entities, these are already
   //- added to the surface
   CubitStatus mesh_surfwoIP(DLIList<TRI *> &facet_list);
   //- mesh the surface with out interior point refinement
   CubitStatus mesh_surfwoIPBC(DLIList<TRI *> &facet_list);
   //- mesh the surface with out interior point refinement or boundary constraints


private:
  ParamTool *pTool; //paramtool needed for sizing function interpolation
  SURF *refFacePtr; // id for surface
  DLIList<TRI *> *facetList; // list of tris created
  DLIList<NODE *> *boundingNodes;
  NODE **boundaryEdgeStartNodes;
  NODE **boundaryEdgeEndNodes; 
  SIZEFUNC *sizingFunction;
  int numBoundaryEdges;
   //- bounding edges on the surfaces represented by their start
   //- and end nodes

  NODE *boxNodes[4];
   //- temporary nodes created enclosing the surface to facilitate
   //- Delaunay node insertion

  int curVisitFlag;
   //- current visited flag for triangles

  CubitStatus insert_hard_points(DLIList<NODE*> &hard_points);
   //- insert hard points into the mesh
  
  CubitStatus init_box();
   //- create two initial bounding triangles

  CubitStatus insert_nodes( DLIList<NODE *> *&bounding_nodes );
   //- insert a list of nodes into an existing triangulation

  CubitStatus constrain_boundary(void);
   //- recover the boundary edges from the mesh

  CubitStatus delete_exterior(void);
   //- delete the triangles on the outside of the boundary and on
   //- the interior of any holes

  CubitStatus delete_tris(DLIList<TRI *> &tri_list, int marked_flag);
   //- delete the triangles in the list with the spacified marked_flag

  CubitStatus mark_tris(TRI *tri_ptr, int &num_found_edges);
   //- recursive function to mark all the triangles up to a marked edge

  CubitStatus refine_interior(DLIList<TRI*> &tri_list);
   //- generate nodes on the interior of the tiangulation to
   //- improve quality. This inserts nodes at the circumcenters of 
   //- triangles roughly based on the algorithm by Jonathon Shewchuk

  CubitStatus classify_triangles(DLIList<TRI *> &tri_list,
                                 DLIList<TRI*>* sorted_lists,
                                 const int num_lists,
                                 const double interval);
   //- order the triangles in the current mesh by their worst angle

  TRI *next_triangle(DLIList<TRI*>* sorted_lists,
                     const int num_lists);
   //- get the next triangle to process;
  
  void clean_up_data(void);
  //- clean up any data we've allocated
  
};

#include "FacetorTool.cpp"
#endif // FACETOR_TOOL_HPP

