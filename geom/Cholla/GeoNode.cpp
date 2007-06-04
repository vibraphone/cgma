//-------------------------------------------------------------------------
// Filename      : GeoNode.cpp
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


#include "GeoNode.hpp"
#include "GeoTet.hpp"


//--------------------------------------------------------------------------
// Function: GeoNode
// Description: constructor
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
GeoNode::GeoNode( CubitVector &xx ) : ownerPtr(NULL)
{
  mLocation = xx;
  isMarked = 0;
}

//--------------------------------------------------------------------------
// Function: GeoNode
// Description: destructor
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
GeoNode::~GeoNode()
{
}

//--------------------------------------------------------------------------
// Function: edge_between
// Description: returns whether an edge exists between the two GeoNodes
// Author: sjowen
// Date: 1/26/2004
//---------------------------------------------------------------------------
CubitBoolean GeoNode::edge_between( GeoNode *other_node )
{
  int ii;
  GeoTet *gtet_ptr;
  for (ii=0; ii<mTetList.size(); ii++)
  {
    gtet_ptr = mTetList.get_and_step();
    if (gtet_ptr->node_index( other_node ) >= 0)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//--------------------------------------------------------------------------
// Function: face_between
// Description: returns whether a face exists between the three GeoNodes
// Author: sjowen
// Date: 1/27/2004
//---------------------------------------------------------------------------
CubitBoolean GeoNode::face_between( GeoNode *other_node0, GeoNode *other_node1 )
{
  int ii;
  GeoTet *gtet_ptr;
  for (ii=0; ii<mTetList.size(); ii++)
  {
    gtet_ptr = mTetList.get_and_step();
    if (gtet_ptr->node_index( other_node0 ) >= 0 &&
        gtet_ptr->node_index( other_node1 ) >= 0)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//--------------------------------------------------------------------------
// Function: tets_at_edge
// Description: get a list of tets shared by the two nodes
//              returns number of tets in list
// Author: sjowen
// Date: 1/27/2004
//---------------------------------------------------------------------------
int GeoNode::tets_at_edge( GeoNode *other_node,                 
                           DLIList<GeoTet *> &gtet_list )
{
  int ntets = 0;
  int ii;
  GeoTet *gtet_ptr;
  for (ii=0; ii<mTetList.size(); ii++)
  {
    gtet_ptr = mTetList.get_and_step();
    if (gtet_ptr->node_index( other_node ) >= 0)
    {
      gtet_list.append( gtet_ptr );
      ntets++;
    }
  }
  return ntets;
}

//--------------------------------------------------------------------------
// Function: interior_edges
// Description: get a list of opposite geonodes to this node
//              don't return nodes that are external to the domain
// Author: sjowen
// Date: 3/14/2004
//---------------------------------------------------------------------------
int GeoNode::interior_edges( DLIList <GeoNode *> &node_list )
{
  int nedges = 0;
  int ii, jj;
  GeoNode *node[4];
  GeoTet *gtet_ptr;
  for (ii=0; ii<mTetList.size(); ii++)
  {
    gtet_ptr = mTetList.get_and_step();
    if (gtet_ptr->inside())
    {
      gtet_ptr->tet_nodes(node[0], node[1], node[2], node[3]);
      for(jj=0; jj<4; jj++)
        node[jj]->marked(1);
    }
  }
  isMarked = 0;

  for (ii=0; ii<mTetList.size(); ii++)
  {
    gtet_ptr = mTetList.get_and_step();
    if (gtet_ptr->inside())
    {
      gtet_ptr->tet_nodes(node[0], node[1], node[2], node[3]);
      for(jj=0; jj<4; jj++)
      {
        if (node[jj]->marked())
        {
          node[jj]->marked(0);
          node_list.append(node[jj]);
          nedges++;
        }
      }
    }
  }
  return nedges; 
}


// EOF

