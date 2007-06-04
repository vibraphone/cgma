//-------------------------------------------------------------------------
// Filename      : GeoTet.hpp
//
// Purpose       : Tet class used for geometry operations.  See also GeoNode.
//
// Description   : light-weight entity used for computational geometry 
//                 tools.  Don't load this up with extra stuff!
//
// Creator       : Steve Owen
//
// Creation Date : 8/13/2003
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------


#ifndef GEOTET_HPP
#define GETTET_HPP

#include "MemoryManager.hpp"
#include "ToolDataUser.hpp"
#include "CubitVector.hpp"

class GeoNode;

class GeoTet : public ToolDataUser
{
private:

  GeoNode *mNodes[4];
    // array ordered nodes on this tet

  cBit isMarked : 1;
    // generic flag for marking the tet

  cBit isInside : 1;
    // defines whether this tet is inside the domain

public:

  GeoTet( GeoNode *nodes[4] );
    //constructor

  ~GeoTet();
    //destructor

  void tet_nodes( GeoNode *&na, GeoNode *&nb, GeoNode *&nc, GeoNode *&nd );
    // get the nodes from this tet

  void tet_face_nodes ( int face_index, GeoNode *&na, GeoNode *&nb, GeoNode *&nc );
    // get the nodes on a specific face of this tet.

  GeoTet *get_connected_tet( GeoNode *na, GeoNode *nb, GeoNode *nc );

  GeoTet *get_connected_tet ( int face_index );

    // get the adjacent tet

  int node_index( GeoNode *node_ptr );
    // return the index of the node in the tet

  void opposite_edge( GeoNode *a_node, GeoNode *b_node, 
                      GeoNode *&c_node, GeoNode *&d_node );
    // return nodes on the opposite edge of the tet.  a_node and 
    // b_node must be on this tet.  c_node and d_node are in
    // no particular order

  void marked( int flag ) { isMarked = flag; }
  int marked() { return isMarked; }
    // get and set the generic marked flag

  void inside( int is_inside ) { isInside = is_inside; }
  int inside() { return isInside; }
    // get and set the isInside flag that defines whether the tet is inside the domain

  CubitStatus circumsphere( CubitVector &circumcenter, double *radius = NULL );
};

// debug functions
void dgtet( GeoTet *tet );
void dgnode( GeoNode *node );

#endif


// EOF

