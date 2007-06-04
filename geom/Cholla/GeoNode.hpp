//-------------------------------------------------------------------------
// Filename      : GeoNode.hpp
//
// Purpose       : Node class used for geometry operations.  See also GeoTet.
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

#ifndef GEONODE_HPP
#define GEONODE_HPP

#include "MemoryManager.hpp"
#include "ToolDataUser.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"

class GeoTet;

class GeoNode: public ToolDataUser
{
private:

  CubitVector mLocation;
    // the location of the node

  DLIList<GeoTet *> mTetList;
    // adjacent tets to this node

  void *ownerPtr;
    // generic owner pointer

  int isMarked;

  static MemoryManager memoryManager;
    //- memory management object

public:

  GeoNode( CubitVector &xx );
    //- constructor

  ~GeoNode();
    //- destructor

  CubitVector coordinates()
    {return mLocation;}
    //- get location of this node

  void set_coordinates( CubitVector &coords )
    { mLocation = coords; }
    //- set coordinate location

  int number_tets()
    {return mTetList.size();}
    //- return number of tets attached to this node.

  void add_tet( GeoTet *tet_ptr, int = CUBIT_TRUE)
    { mTetList.append(tet_ptr); }
    //- add a tet to the adjacency list

  void remove_tet( GeoTet *tet_ptr )
    { mTetList.remove(tet_ptr); }
    //- remove tet from the adjacency list

  DLIList<GeoTet*> *tet_list()
    { return &mTetList; }
    //- get the pointer to the tet list

  void *get_owner() {return ownerPtr;}
  void set_owner( void *owner_ptr ) {ownerPtr = owner_ptr;}
    //- get and set the generic owner pointer

  void marked( int flag ) { isMarked = flag; }
  int marked() { return isMarked; }
    //- get and set the marked flag

  CubitBoolean edge_between( GeoNode *other_node );
    //- returns whether an edge exists between the two GeoNodes

  CubitBoolean face_between( GeoNode *other_node0, GeoNode *other_node1 );
    //- returns whether a face exists between the three GeoNodes

  int tets_at_edge( GeoNode *other_node,                 
                    DLIList<GeoTet *> &gtet_list );
    //- get a list of tets shared by the two nodes
    //- returns number of tets in list

  int interior_edges( DLIList <GeoNode *> &node_list );
    //- get a list of opposite geonodes to this node
    //- don't return nodes that are external to the domain

};

#endif

// EOF

