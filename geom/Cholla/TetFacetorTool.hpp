//-------------------------------------------------------------------------
// Class         : TetFacetorTool
//
// Filename      : TetFacetorTool.hpp
//
// Purpose       : 3D Delaunay Tesselator.  Given a set of 3D points, will define 
//                 a Delaunay tesselation.  Does not currently support boundary 
//                 constraints.  Result will be the convex hull of the points
//
// Description   : TetFacetorTool is a template class that currently accepts 
//                 MeshEntity classes.  Any other data structure that satisfies
//                 the appropriate member functions should also be handled.
//                 see also FacetorTool  
//
//                 Following is the list of required members of the TET and 
//                 NODE classes to use this template
//                 
//                 TET: add_TD, get_TD, tet_nodes, get_connected_tet, 
//                      new, delete
//                 NODE: number_edges, coordinates, add_tet, remove_tet
//                       new, delete
//
// Creator       : Steve Owen
//
// Creation Date : 8/2/2003
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

#ifndef TET_FACETOR_TOOL_HPP
#define TET_FACETOR_TOOL_HPP

//#define DEBUG_TET_FACETOR

#include "CubitBox.hpp"

#ifdef DEBUG_TET_FACETOR
#define TET CubitTet
#define SUBTET NodeTet
#define NODE CubitNode
class TET;
class SUBTET;
class NODE;
#else
template<class TET, class SUBTET, class NODE>
#endif 
class TetFacetorTool
{
public:
  TetFacetorTool();
    // constructor

  ~TetFacetorTool(void);
    // destructor

  int initialize(CubitBox &bounding_box);
    //- initialize the Delaunay domain with initial tets.  Use if
    //- point-by-point insertion will be performed.  (Don't use
    //- with tesselate())

  int insert_one_node( NODE *node_ptr );
    //- insert a single node into tesselation

  int finish();
    //- complete the tesselation - remove bounding box nodes

  int get_tets( std::vector<TET *> &tet_list );
    //- retreive the current list of tets 

  TET *get_outside_tet();
    //- get one tet that is on the outside of the boundary

  int get_interior_tets( std::vector<TET *> &tet_list );
    //- get tets that are not connected to a bounding box node

  int tesselate(std::vector<NODE *> &node_list, std::vector<TET *> &tet_list);
    //- given a list of points, form a Delaunay tesselation.  Does
    //- initializion and finish (All-in-one function)

  int read_data( const char *filename, std::vector<NODE *>&node_list );
    //- read data from a file and create list of NODEs

  int circumsphere( TET *tet_ptr, CubitVector &center, double &radius2 );
    // get the circumsphere info for the tet

  int natural_neighbor_tets( CubitVector &xx, std::vector<TET *> &neighbor_tet_list, 
                             NODE *&duplicate_node );
    // get all tets whose circumsphere contain the point

  int watson_insert( NODE *node_ptr, std::vector<TET *> &neighbor_tet_list );
    // insert using Bowyer-Watson algrithm

  double get_tol(){return csTol;}
  void set_tol(double tol){csTol=tol;}
    // get and set circumsphere comparison tolerance.  csTol is 
    // computed by default in create_bbox_tets

private:

  std::vector<TET *> tetList;
    // current set of tets in the set

  TET *lastTet;
    // last tet inserted/visited

  int tetVisited;
    // flag keeps track current iteration of Delaunay insertion

  CubitBox bBox;
    // the bounding box of the nodes

  double csTol;
    // circumsphere tolerance

  NODE *boxNodes[8];
    // nodes on the bounding box of the tet

  int mDebug;
    // debug flag;

  int init_box(std::vector<NODE *> &node_list);
    // create the initial tets where all subsequent nodes will be inserted

  int create_bbox_tets();
    // create the initial tets for the Delaunay domain.

  int remove_bbox_tets();
    // remove the tets at the bounding box

  int is_bbox(NODE *n);
    // return if this is a node on the bounding box;

  int insert_nodes( std::vector<NODE *> &node_list );
    // insert nodes into existing tesselation

  int locate_point( CubitVector &xx, NODE *&exact_node, TET *&containing_tet );
    // return the tet the point is located in.  Do a walking
    // algorithm starting from the lastTet.

  int exhaustive_locate_point( CubitVector &xx, NODE *&exact_node, 
                               TET *&containing_tet );
    // Try all tets in the list (that haven't already been tried)
    // This is called only when locate_point fails

  CubitBoolean point_in_circumsphere( TET *adj_tet, CubitVector &xx, 
                                     std::vector<TET *> &neighbor_tet_list );
    // recursive function, determines if the point is within the
    // circumsphere of the given tet and then recurses to its
    // neighbors if it is.  Add to the neighbor_tet_list as we go.

  void set_tet_visited( TET *tet_ptr, int new_visit_flag );
    // set the tet visited flag

  int tet_visited( TET *tet_ptr  );
    // return the visited flag

  double tet_volume( const CubitVector &a, 
                     const CubitVector &b, 
                     const CubitVector &c, 
                     const CubitVector &d );
    // compute signed volume of tet
  
};

#ifndef DEBUG_TET_FACETOR
#include "TetFacetorTool.cpp"
#endif


#endif // TET_FACETOR_TOOL_HPP

