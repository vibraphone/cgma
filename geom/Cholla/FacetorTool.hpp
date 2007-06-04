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

  DLIList<TRI*> *triSortArray;
   // array of lists of triangles sorted by their worst angle

  NODE *boxNodes[4];
   //- temporary nodes created enclosing the surface to facilitate
   //- Delaunay node insertion
    
  TRI *lastTriLocated;
  //TRI *debugTri;
   //- last triangle the locate point function found

  int curVisitFlag;
   //- current visited flag for triangles

  CubitStatus insert_hard_points(DLIList<NODE*> &hard_points);
   //- insert hard points into the mesh
  
  CubitStatus init_box();
   //- create two initial bounding triangles

  CubitVector &circumcenter( TRI *tri_ptr );
   //- get the circumcenter of the triangle

  double radius( TRI *tri_ptr );
   //- get the radius squared of the triangle circumcircle

  CubitBoolean tri_visited( TRI *tri_ptr );
  void tri_visited( TRI *tri_ptr, CubitBoolean visited );
   //- get and set the visited flag

  void tri_sort_list( TRI *tri_ptr, int sort_list_index );
  int tri_sort_list( TRI *tri_ptr );
   //- get and set the tri sort array index that this tri is in

  CubitStatus insert_nodes( DLIList<NODE *> *&bounding_nodes );
   //- insert a list of nodes into an existing triangulation

  CubitStatus insert_node( NODE *node_ptr );
    //- nsert one node into Delaunay mesh          

  CubitBoolean are_nodes_colinear( TRI *tri_ptr );
   //- are the nodes of this tri roughly colinear
  
  CubitStatus bowyer_watson_insert( NODE *node_ptr,
                                    DLIList<TRI *> &tri_list);
   //- Bowyer-Watson insertion into an existing Delaunay Mesh
  CubitStatus valid_void( NODE *node_ptr,
                          DLIList<TRI *> &tri_list);
    //check the void created in Bowyer-Watson to ensure the boundary
    //seems valid.

  CubitStatus natural_neighbor_tris( CubitVector &the_point,
                                     DLIList <TRI *> &tri_list);
   //- get a list of all triangles whose circumcircle contain 
   //- the point

  CubitStatus locate_point( CubitVector &the_point,
                            TRI *&tri_ptr );
   //- return the triangle the point is located in

  CubitStatus exhaustive_locate_point( CubitVector &the_point,
                                       TRI *&tri_ptr );
   //- return the triangle the point is located in by checking 
   //- all triangles

  CubitStatus point_in_circumcircle( CubitVector &the_point,
                                     TRI *tri_ptr,
                                     DLIList <TRI *> &tri_list );
   //- determine if the point is inside the circumcircle of the
   //- triangle and recurse to the adjacent triangles 

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

  CubitStatus classify_triangles(DLIList<TRI *> &tri_list);
   //- order the triangles in the current mesh by their worst angle

  CubitStatus classify_tri_by_angle( TRI *tri_ptr );
   //- compute the angles at the triangle vertices and classify
   //- by its worst triangle

  TRI *next_triangle(void);
   //- get the next triangle to process

  double get_size(CubitVector &cc, TRI *tri_ptr);

  CubitStatus insert_at_circumcenter(TRI *tri_ptr);
   //- insert a node at the circumcenter of a tri
  
  void clean_up_data(void);
  //- clean up any data we've allocated
  
};

// Added by CAT for NT port
#if defined(TEMPLATE_DEFS_INCLUDED)
  #include "FacetorTool.cpp"
#endif
#endif // FACETOR_TOOL_HPP

