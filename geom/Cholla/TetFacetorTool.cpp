//-------------------------------------------------------------------------
// Filename      : TetFacetorTool.cpp
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
//                 NODE: number_tets, coordinates, add_tet, remove_tet
//                       new, delete, tet_list
//
// Creator       : Steve Owen
//
// Creation Date : 8/2/2003
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

//#define DEBUG_TET_FACETOR


#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

// defining DEBUG_TET_FACETOR will allow you to compile and link as a 
// standard C++ class rather than a template.  The definitions used for
// NODE and TET are the Cubit mesh entities.  Since they can't be used
// in the geom directory, DEBUG_TET_FACETOR should not be defined under
// normal circumstances.
#ifdef DEBUG_TET_FACETOR
#include "..\mdb\CubitTet.hpp"
#include "..\mdb\NodeTet.hpp"
#include "..\mdb\CubitNode.hpp" 
#define TET CubitTet
#define SUBTET NodeTet
#define NODE CubitNode
#define TET_FACETOR_TOOL    TetFacetorTool::
#define TET_FACETOR_TOOL__I int TetFacetorTool::
#define TET_FACETOR_TOOL__V void TetFacetorTool::
#define TET_FACETOR_TOOL__D double TetFacetorTool::
#define TET_FACETOR_TOOL__B CubitBoolean TetFacetorTool::
#define TET_FACETOR_TOOL__T TET *TetFacetorTool::
#else
#define TET_FACETOR_TOOL    template<class TET, class SUBTET, class NODE> MY_INLINE TetFacetorTool<TET, SUBTET, NODE>::
#define TET_FACETOR_TOOL__I template<class TET, class SUBTET, class NODE> MY_INLINE int TetFacetorTool<TET, SUBTET, NODE>::
#define TET_FACETOR_TOOL__V template<class TET, class SUBTET, class NODE> MY_INLINE void TetFacetorTool<TET, SUBTET, NODE>::
#define TET_FACETOR_TOOL__D template<class TET, class SUBTET, class NODE> MY_INLINE double TetFacetorTool<TET, SUBTET, NODE>::
#define TET_FACETOR_TOOL__B template<class TET, class SUBTET, class NODE> MY_INLINE CubitBoolean TetFacetorTool<TET, SUBTET, NODE>::
#define TET_FACETOR_TOOL__T template<class TET, class SUBTET, class NODE> MY_INLINE TET *TetFacetorTool<TET, SUBTET, NODE>::
#endif

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif
#include <math.h>
#include "GeometryDefines.h"
#include "TetFacetorTool.hpp"
#include "TDDelaunay.hpp"
#include "TDInterpNode.hpp"


//-------------------------------------------------------------------------
// Function:    TetFacetorTool
// Description: constructor
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL 
TetFacetorTool(  )
{
  //update private variables
  mDebug = 1;

}


//-------------------------------------------------------------------------
// Function:    ~TetFacetorTool
// Description: destructor
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL
~TetFacetorTool()
{
}

//-------------------------------------------------------------------------
// Function:    initialize
// Description: initialize the Delaunay domain based on the bounding box
//              provided.  All nodes to be inserted later are intended to
//              lie within this box.  Note that the bounding box should
//              be larger than bounding box of nodes by at least about 10-20%
// Author:      sjowen
// Date:        1/22/2004
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
initialize(CubitBox &bounding_box)
{
  bBox = bounding_box;
  return create_bbox_tets();
}

//-------------------------------------------------------------------------
// Function:    finish
// Description: complete the tesselation.  Remove bounding box nodes and
//              their connected tets
// Author:      sjowen
// Date:        1/22/2004
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
finish()
{
  return remove_bbox_tets();
}

//-------------------------------------------------------------------------
// Function:    get_tets
// Description: get the current tet list
// Author:      sjowen
// Date:        1/22/2004
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
get_tets(DLIList<TET*> &tet_list)
{
  tet_list += tetList;
  return tet_list.size();
}

//-------------------------------------------------------------------------
// Function:    get_outside_tet
// Description: get the current tet list
// Author:      sjowen
// Date:        1/22/2004
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__T
get_outside_tet()
{
  NODE *node = boxNodes[0];
  DLIList<TET *> *tet_list_ptr = node->tet_list();
  assert(tet_list_ptr != NULL && tet_list_ptr->size() > 0);
  tet_list_ptr->reset();
  return tet_list_ptr->get();
}

//-------------------------------------------------------------------------
// Function:    get_interior_tets
// Description: get tets that are not connected to a bounding box node
// Author:      sjowen
// Date:        1/22/2004
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
get_interior_tets(DLIList<TET*> &tet_list)
{
  int ii, jj;
  TET *tet_ptr;
  NODE *nodes[4];
  int found = 0;
  for (ii=0; ii<tetList.size(); ii++)
  {
    tet_ptr = tetList.get_and_step();  
    tet_ptr->tet_nodes(nodes[0], nodes[1], nodes[2], nodes[3]);
    found = 0;
    for(jj=0; jj<4 && !found; jj++)
    {
      if (is_bbox(nodes[jj]))
      {
        found = 1;
      }
    }
    if (!found)
      tet_list.append(tet_ptr);
  }

  return tet_list.size();
}

//-------------------------------------------------------------------------
// Function:    tesselate
// Description: given a set of points, create delaunay tets from them
//              does all initialization and removal of bounding box tets
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
tesselate(DLIList<NODE *> &node_list, DLIList<TET *> &tet_list)
{
  int stat;

  // create a set of tets in a bounding box to start
  stat = init_box( node_list );

  // insert all the nodes
  if (stat == CUBIT_SUCCESS)
  {
    int num_failed = insert_nodes( node_list );
    stat = (num_failed==0) ? CUBIT_SUCCESS : CUBIT_FAILURE;
  }

  // remove the nodes at the bounding box.
  // this should give you a convex hull of the points (provided your points
  // are dense enough)

  if (stat == CUBIT_SUCCESS)
  {
    stat = remove_bbox_tets();
  }

  // copy the local list to the tet_list argument to pass back
  if (stat == CUBIT_SUCCESS)
  { 
    tet_list += tetList;
  }

  return stat;
}

//-------------------------------------------------------------------------
// Function:    init_box
// Description: create a set of tets in a bounding box to start
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
init_box(DLIList<NODE *> &node_list)
{
  CubitVector minbox(CUBIT_DBL_MAX, CUBIT_DBL_MAX, CUBIT_DBL_MAX);
  CubitVector maxbox(-CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX);

  // loop through nodes to set bounding box

  CubitVector coor;
  NODE *node_ptr;
  int ii;
  for (ii=0; ii<node_list.size(); ii++)
  {
    node_ptr = node_list.get_and_step();
    coor = node_ptr->coordinates();

    if (coor.x() > maxbox.x()) maxbox.x( coor.x() );
    if (coor.x() < minbox.x()) minbox.x( coor.x() );
    if (coor.y() > maxbox.y()) maxbox.y( coor.y() );
    if (coor.y() < minbox.y()) minbox.y( coor.y() );
    if (coor.z() > maxbox.z()) maxbox.z( coor.z() );
    if (coor.z() < minbox.z()) minbox.z( coor.z() );
  }

  // Expand the bounding box by 20% of the size of the diagonal

  double dx = maxbox.x() - minbox.x();
  double dy = maxbox.y() - minbox.y();
  double dz = maxbox.z() - minbox.z();
  double expand = 0.2 * sqrt(SQR(dx) + SQR(dy) + SQR(dz));
  minbox.x( minbox.x() - expand );
  minbox.y( minbox.y() - expand );
  minbox.z( minbox.z() - expand );
  maxbox.x( maxbox.x() + expand );
  maxbox.y( maxbox.y() + expand );
  maxbox.z( maxbox.z() + expand );

  bBox.reset( minbox, maxbox );

  return create_bbox_tets();
}

//-------------------------------------------------------------------------
// Function:    create_bbox_tets
// Description: create the tets that contain a node on the bounding box
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
create_bbox_tets()
{
  // use the class global bBox as the bounds of the initial tets

  CubitVector maxbox, minbox;
  maxbox = bBox.maximum();
  minbox = bBox.minimum();

  // Initialize the class global vars

  lastTet = NULL;
  tetVisited = CUBIT_INT_MIN;

  // For tolerance, find a representative (non-zero) number from
  // the nodes to determine relative magnitude of numbers.  Take the
  // log of this number, subtract 6 from it, then use this number
  // as the exponent (likely a negative number) for the tolerance

  double tol = CUBIT_MAX(fabs(maxbox.x()),fabs(maxbox.y()));
  tol = CUBIT_MAX(tol,fabs(maxbox.z()));
  double temp = CUBIT_MAX(fabs(minbox.x()),fabs(minbox.y()));
  temp = CUBIT_MAX(temp,fabs(minbox.z()));
  tol = CUBIT_MAX(tol,temp);
  tol = pow(10.0, (double)(log10(tol) - 6.0));
  csTol = tol * tol;

  // create the bounding box nodes

  CubitVector coor;
  int ii;
  for (ii=0; ii<8; ii++) 
  {
    switch (ii) {
      case 0: coor.set(minbox.x(), minbox.y(), minbox.z()); break;
      case 1: coor.set(maxbox.x(), minbox.y(), minbox.z()); break;
      case 2: coor.set(maxbox.x(), maxbox.y(), minbox.z()); break;
      case 3: coor.set(minbox.x(), maxbox.y(), minbox.z()); break;
      case 4: coor.set(minbox.x(), minbox.y(), maxbox.z()); break;
      case 5: coor.set(maxbox.x(), minbox.y(), maxbox.z()); break;
      case 6: coor.set(maxbox.x(), maxbox.y(), maxbox.z()); break;
      case 7: coor.set(minbox.x(), maxbox.y(), maxbox.z()); break;
    }
    boxNodes[ii] = new NODE(coor);
  }

  // create 5 tets to fill the box
  
  const int itet_config[5][4] = {{0,1,3,4},{1,2,3,6},{1,5,6,4},
                                 {3,7,4,6},{1,6,3,4}};
  NODE *tet_nodes[4];
  TET *new_tet;
  int jj;
  for (ii=0; ii<5; ii++)
  {
    for (jj=0; jj<4; jj++)
    {
      tet_nodes[jj] = boxNodes[itet_config[ii][jj]];
    }
    new_tet = new SUBTET( tet_nodes );
    for (jj=0; jj<4; jj++)
      tet_nodes[jj]->add_tet(new_tet);
    tetList.append( new_tet );
  }


  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    remove_bbox_tets
// Description: remove the tets that contain a node on the bounding box
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
remove_bbox_tets()
{
  int stat = CUBIT_SUCCESS;
  int found;
  TET *tet_ptr;
  int ii, jj, kk;

  NODE *nodes[4];
  for (ii=0; ii<tetList.size(); ii++)
  {
    found = 0;
    tet_ptr = tetList.get();  
    tet_ptr->tet_nodes(nodes[0], nodes[1], nodes[2], nodes[3]);
    for(jj=0; jj<4 && !found; jj++)
    {
      if (is_bbox(nodes[jj]))
      {
        found = 1;
        for(kk=0; kk<4; kk++)
          nodes[kk]->remove_tet(tet_ptr);
        tetList.change_to(NULL);
        delete tet_ptr;
      }
    }
    tetList.step();
  }
  tetList.remove_all_with_value(NULL);

  for (ii=0; ii<8; ii++)
  {
    delete boxNodes[ii];
  }

  return stat;
}

//-------------------------------------------------------------------------
// Function:    is_bbox
// Description: determine if this node is one of the bounding box nodes
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
is_bbox(NODE *n)
{
  int ii;
  for(ii=0; ii<8; ii++)
    if(n == boxNodes[ii])
      return 1;
  return 0;
}

//-------------------------------------------------------------------------
// Function:    insert_nodes
// Description: insert nodes into existing tesselation (aasumes at least
//              bounding box tets exist)
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
insert_nodes( DLIList<NODE *> &node_list )
{
  int ii;
  NODE *node_ptr;
  int num_failed = 0;
  for (ii=0; ii<node_list.size(); ii++)
  {
    node_ptr = node_list.get_and_step();
    if(insert_one_node( node_ptr ) == CUBIT_FAILURE)
    {
      num_failed++;
    }
  }

  // If any failed to insert the first time, then try them again

  if (num_failed > 0) 
  {
    num_failed = 0;
    for (ii=0; ii<node_list.size(); ii++)
    {
      node_ptr = node_list.get_and_step();
      if (node_ptr->number_tets() == 0)
      {
        if(insert_one_node( node_ptr ) == CUBIT_FAILURE)
        {
         num_failed++;
        }
      }
    }
  }
      
  return num_failed;

}

//-------------------------------------------------------------------------
// Function:    insert_one_node
// Description: insert a single node into the tesselation. (Must be within 
//              existing convex hull)
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
insert_one_node( NODE *node_ptr )
{

  // Get all tets whose circumsphere contain the point

  CubitVector xx = node_ptr->coordinates();
  DLIList<TET *> neighbor_tet_list;
  NODE *duplicate_node = NULL;
  int stat = natural_neighbor_tets( xx, neighbor_tet_list, duplicate_node );
  if (stat != CUBIT_SUCCESS)
    return stat;


  // If this new node fell exactly on top of an existing node, then
  // ignore it

  if (duplicate_node != NULL)
  {
    PRINT_WARNING("Duplicate node detected in Delaunay insertion at (%f %f %f).\n"
                  "  Ignoring node and continuing.\n", xx.x(), xx.y(), xx.z());
    return CUBIT_SUCCESS;
  }

  // insert the node using the bowyer-watson algorithm

  stat = watson_insert( node_ptr, neighbor_tet_list );

  return stat;

}

//-------------------------------------------------------------------------
// Function:    natural_neighbor_tets
// Description: Get all tets whose circumsphere contain the point
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
natural_neighbor_tets( CubitVector &xx, DLIList<TET *> &neighbor_tet_list, 
                       NODE *&exact_node )
{
  // Determine the tet where the point is located.  If its at a node
  // return now with success (trivial interpolation case)

  TET *containing_tet;
  int stat = locate_point( xx, exact_node, containing_tet );
  if (stat != CUBIT_SUCCESS)
    return stat;
  if (exact_node != NULL)
    return CUBIT_SUCCESS;

  // Put the tet that contains the point as the first one on the list 
  // and mark it as visited

  neighbor_tet_list.append( containing_tet );
  set_tet_visited( containing_tet, ++tetVisited );

  // Recursively search, (starting with the tet the point is in)
  // search for all tets whose circumsphere contain the point and place 
  // on the neighbor_tet_list

  int iface;
  TET *adj_tet = NULL;
  NODE *na, *nb, *nc, *nd;
  containing_tet->tet_nodes(na, nb, nc, nd);
  for (iface=0; iface<4; iface++) 
  {
    switch(iface)
    {
      case 0: adj_tet = (TET *)containing_tet->get_connected_tet( nb, nc, nd ); break;
      case 1: adj_tet = (TET *)containing_tet->get_connected_tet( nd, nc, na ); break;
      case 2: adj_tet = (TET *)containing_tet->get_connected_tet( nd, na, nb ); break;
      case 3: adj_tet = (TET *)containing_tet->get_connected_tet( nc, nb, na ); break;
    }
    if (adj_tet != NULL)
    {
      if (tet_visited(adj_tet) != tetVisited)
      {
        point_in_circumsphere( adj_tet, xx, neighbor_tet_list );
      }
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    watson_insert
// Description: insert a single node into the tesselation using the 
//              bowyer watson algorithm
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
watson_insert( NODE *node_ptr, DLIList<TET *> &neighbor_tet_list )
{
  TET *tet_ptr;
  int ii,jj;

  // mark all the tets in the list

  tetVisited++;
  for (ii=0; ii<neighbor_tet_list.size(); ii++)
    set_tet_visited(neighbor_tet_list.get_and_step(), tetVisited);

  // Go through the neighbor tets and find the faces that are not 
  // adjacent to any other tet in the neighbor list.  These faces
  // will serve as base facets for new tets

  DLIList<TET *> new_tet_list;
  CubitVector cc;
  double radsq;
  TET *adj_tet = NULL, *new_tet = NULL;
  DLIList<NODE *> cavity_nodes;
  NODE *nodes[4];
  int itet, iface;
  int ntet = 0;
  NODE *na, *nb, *nc, *nd;

  // form the list of boundary faces of the cavity
  for (itet=0; itet<neighbor_tet_list.size(); itet++) 
  {
    tet_ptr = neighbor_tet_list.get_and_step();
    tet_ptr->tet_nodes(na, nb, nc, nd);
    for (iface=0; iface<4; iface++) 
    {  
      switch(iface)
      {
        case 0: adj_tet = (TET *)tet_ptr->get_connected_tet( nb, nc, nd ); break;
        case 1: adj_tet = (TET *)tet_ptr->get_connected_tet( nd, nc, na ); break;
        case 2: adj_tet = (TET *)tet_ptr->get_connected_tet( nd, na, nb ); break;
        case 3: adj_tet = (TET *)tet_ptr->get_connected_tet( nc, nb, na ); break;
      }
      if (adj_tet == NULL || tet_visited(adj_tet) != tetVisited)
      {
        switch(iface)
        {
          case 0:
            cavity_nodes.append(nb);
            cavity_nodes.append(nd);
            cavity_nodes.append(nc);
            break;
          case 1:
            cavity_nodes.append(na);
            cavity_nodes.append(nc);
            cavity_nodes.append(nd);
            break;
          case 2:
            cavity_nodes.append(na);
            cavity_nodes.append(nd);
            cavity_nodes.append(nb);
            break;
          case 3:
            cavity_nodes.append(na);
            cavity_nodes.append(nb);
            cavity_nodes.append(nc);
            break;
        }
        ntet++;
      }
    }
  }

  cavity_nodes.reset();
  for (itet=0; itet<ntet; itet++)
  {

    // form a new tet with this face and the node
        
    nodes[0] = cavity_nodes.get_and_step();
    nodes[1] = cavity_nodes.get_and_step();
    nodes[2] = cavity_nodes.get_and_step();
    nodes[3] = node_ptr;

    // check volume of this new tet - if its less than zero then the cavity is 
    // not strictly star-shaped.  fail the insertion, delete any new tets already
    // created and return now.

    double vol = tet_volume(nodes[0]->coordinates(), 
                            nodes[2]->coordinates(), 
                            nodes[1]->coordinates(), 
                            nodes[3]->coordinates());
    if (vol < 0.0)
    {
      for (ii=0; ii<new_tet_list.size(); ii++)
      {
        tet_ptr = new_tet_list.get_and_step(); 
        tet_ptr->tet_nodes(nodes[0], nodes[1], nodes[2], nodes[3]);
        for(jj=0; jj<4; jj++)
          nodes[jj]->remove_tet(tet_ptr);
        delete tet_ptr;        
      }
      return CUBIT_FAILURE;
    }
    new_tet = new SUBTET(nodes);
    for (ii=0; ii<4; ii++)
      nodes[ii]->add_tet(new_tet);
    new_tet_list.append(new_tet);

    // define the circumsphere.  If it fails, then fail this
    // insertion, delete any tets already created and return now.

    if (circumsphere(new_tet, cc, radsq) == CUBIT_FAILURE)
    {
      for (ii=0; ii<new_tet_list.size(); ii++)
      {
        tet_ptr = new_tet_list.get_and_step(); 
        tet_ptr->tet_nodes(nodes[0], nodes[1], nodes[2], nodes[3]);
        for(jj=0; jj<4; jj++)
          nodes[jj]->remove_tet(tet_ptr);
        delete tet_ptr;       
      }
      return CUBIT_FAILURE;
    }
    set_tet_visited(new_tet, tetVisited);
  }

  // Set last face for Locate Point

  if (new_tet != NULL)
    lastTet = new_tet;

  // Delete all tets in the neighbor tet list

  for (itet=0; itet<neighbor_tet_list.size(); itet++) 
  {
    tet_ptr = neighbor_tet_list.get_and_step();
    tet_ptr->tet_nodes(nodes[0], nodes[1], nodes[2], nodes[3]);
    for(jj=0; jj<4; jj++)
      nodes[jj]->remove_tet(tet_ptr);
    tetList.move_to_nearby(tet_ptr);
    tetList.extract();
    delete tet_ptr;
  }

  tetList += new_tet_list;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    locate_point
// Description: return the tet the point is located in.  Do a walking
//              algorithm starting from the lastTet.
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I 
locate_point( CubitVector &xx, 
              NODE *&exact_node, 
              TET *&containing_tet )
{
  // Use the last tet to be located as the first try.  If there is
  // no last tet, then use one tet on the global list

  exact_node = NULL;
  containing_tet = NULL;
  TET *cur_tet = lastTet;
  if (cur_tet == NULL)
  {
    cur_tet = tetList.get();
  }

  // keep track of the tets we visit by marking as we go - 
  // increment the visit flag for this call
  tetVisited++;

  TET *adj_tet = NULL;
  double tol = GEOMETRY_RESABS;
  NODE *na, *nb, *nc, *nd;
  CubitVector A, B, C, D;
  double volcoord_A, volcoord_B, volcoord_C, volcoord_D, vol;
  CubitBoolean found = CUBIT_FALSE;
  while (!found) 
  {
    set_tet_visited( cur_tet, tetVisited );

    // Get the coords of the nodes on the tet and compute the
    // barycentric coords (areacoords) of the point with respect 
    // to the tet.  If all area coords are positive, then 
    // the point is inside the tet.  If any are negative
    // then choose the next try as the adjacent tet to
    // the area coord that was the smallest

    cur_tet->tet_nodes( na, nb, nc, nd );
    A = na->coordinates();
    B = nb->coordinates();
    C = nc->coordinates();
    D = nd->coordinates();
    volcoord_A = tet_volume(B, C, D, xx);
    volcoord_B = tet_volume(A, D, C, xx);
    volcoord_C = tet_volume(A, B, D, xx);
    volcoord_D = tet_volume(A, C, B, xx);
    vol = tet_volume(A, C, B, D);

    // check for inverted tet
    assert (vol > 0.0);

    // normalize coords

    //vol = volcoord_A + volcoord_B + volcoord_C + volcoord_D;
    if (fabs(vol) > CUBIT_RESABS) 
    {
      volcoord_A /= vol;
      volcoord_B /= vol;
      volcoord_C /= vol;
      volcoord_D /= vol;
      tol = GEOMETRY_RESABS;
    }
    else 
    {
      tol = csTol;
    }

    if (volcoord_A > -tol && 
        volcoord_B > -tol && 
        volcoord_C > -tol &&
        volcoord_D > -tol) 
    {
      found = TRUE;

      // if three of the areacoords are +-tol then we are at an existing node

      if (volcoord_B < tol && volcoord_C < tol && volcoord_D < tol) 
      {
        exact_node = na;
      }
      else if (volcoord_A < tol && volcoord_C < tol && volcoord_D < tol) 
      {
        exact_node = nb;
      }
      else if (volcoord_A < tol && volcoord_B < tol && volcoord_D < tol) 
      {
        exact_node = nc;
      }
      else if (volcoord_A < tol && volcoord_B < tol && volcoord_C < tol) 
      {
        exact_node = nd;
      }

      // this is the general case where all areacoords are positive.
      // we have located the tet xx is contained in

      else
      {
        containing_tet = cur_tet;
      }
    }

    // at least one of the areacoords is negative.  Advance to the next adjacent
    // tet.  Choose the adjacent tet based on the sign of the areacoords

    else 
    {
      if (volcoord_A <= volcoord_B && 
          volcoord_A <= volcoord_C &&
          volcoord_A <= volcoord_D) 
      {
        adj_tet = (TET *)cur_tet->get_connected_tet( nb, nc, nd );
      }
      else if (volcoord_B <= volcoord_A && 
               volcoord_B <= volcoord_C &&
               volcoord_B <= volcoord_D) 
      {
        adj_tet = (TET *)cur_tet->get_connected_tet( nd, nc, na );
      }
      else if (volcoord_C <= volcoord_A && 
               volcoord_C <= volcoord_B &&
               volcoord_C <= volcoord_D) 
      {
        adj_tet = (TET *)cur_tet->get_connected_tet( nd, na, nb );
      }
      else 
      {
        adj_tet = (TET *)cur_tet->get_connected_tet( nc, nb, na );
      }
      cur_tet = adj_tet;
      
      // Check if we just left the triangulation or we have already been 
      // to this tet.  If so, use the exhaustive search
      
      if(cur_tet == NULL  || tet_visited( cur_tet ) == tetVisited)
      {
        return exhaustive_locate_point(xx, exact_node, containing_tet);
      }
    }
  }
  lastTet = containing_tet;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    exhaustive_locate_point
// Description: Try all tets in the list (that haven't already been tried)
//              This is called only when locate_point fails
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I 
exhaustive_locate_point( CubitVector &xx, 
              NODE *&exact_node, 
              TET *&containing_tet )
{
  // Use the last tet to be located as the first try.  If there is
  // no last tet, then use one tet on the global list

  exact_node = NULL;
  containing_tet = NULL;
  TET *cur_tet;
  double tol = GEOMETRY_RESABS;
  NODE *na, *nb, *nc, *nd;
  CubitVector A, B, C, D;
  double volcoord_A, volcoord_B, volcoord_C, volcoord_D, vol;
  CubitBoolean found = CUBIT_FALSE;
  for (int ii=0; ii<tetList.size() && !found; ii++)
  {
    cur_tet = tetList.get_and_step();

    // Don't chak tets we have already chaked.  tetVisited should be 
    // current (was set in locate_point)

    if (tet_visited(cur_tet) == tetVisited)
      continue;
    set_tet_visited( cur_tet, tetVisited );

    // Get the coords of the nodes on the tet and compute the
    // barycentric coords (areacoords) of the point with respect 
    // to the tet.  If all area coords are positive, then 
    // the point is inside the tet.  If any are negative
    // then choose the next try as the adjacent tet to
    // the area coord that was the smallest

    cur_tet->tet_nodes( na, nb, nc, nd );
    A = na->coordinates();
    B = nb->coordinates();
    C = nc->coordinates();
    D = nd->coordinates();
    volcoord_A = tet_volume(B, C, D, xx);
    volcoord_B = tet_volume(A, D, C, xx);
    volcoord_C = tet_volume(A, B, D, xx);
    volcoord_D = tet_volume(A, C, B, xx);
    vol = tet_volume(A, C, B, D);

    // normalize coords

    //vol = volcoord_A + volcoord_B + volcoord_C + volcoord_D;
    if (fabs(vol) > CUBIT_RESABS) 
    {
      volcoord_A /= vol;
      volcoord_B /= vol;
      volcoord_C /= vol;
      volcoord_D /= vol;
      tol = GEOMETRY_RESABS;
    }
    else 
    {
      tol = csTol;
    }

    if (volcoord_A > -tol && 
        volcoord_B > -tol && 
        volcoord_C > -tol &&
        volcoord_D > -tol) 
    {
      found = TRUE;

      // if three of the areacoords are +-tol then we are at an existing node

      if (volcoord_B < tol && volcoord_C < tol && volcoord_D < tol) 
      {
        exact_node = na;
      }
      else if (volcoord_A < tol && volcoord_C < tol && volcoord_D < tol) 
      {
        exact_node = nb;
      }
      else if (volcoord_A < tol && volcoord_B < tol && volcoord_D < tol) 
      {
        exact_node = nc;
      }
      else if (volcoord_A < tol && volcoord_B < tol && volcoord_C < tol) 
      {
        exact_node = nd;
      }

      // this is the general case where all areacoords are positive.
      // we have located the tet xx is contained in

      else
      {
        containing_tet = cur_tet;
      }
    }
  }
  if (!found)
  {
    lastTet = NULL;
    return CUBIT_FALSE;
  }
  lastTet = containing_tet;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    point_in_circumsphere
// Description: recursive function, determines if the point is within the
//              circumsphere of the given tet and then recurses to its
//              neighbors if it is.  Add to the neighbor_tet_list as we go.
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__B 
point_in_circumsphere( TET *tet_ptr, 
                       CubitVector &xx, 
                       DLIList<TET *> &neighbor_tet_list )
{
  // The circumsphere of the face has been previously stored with the
  // face.  Determine distance squared from center of circle to point

  set_tet_visited(tet_ptr, tetVisited);
  double radsq;
  CubitVector cc;
  circumsphere(tet_ptr, cc, radsq);

  double dist2 = (SQR(xx.x()-cc.x()) + 
                  SQR(xx.y()-cc.y()) + 
                  SQR(xx.z()-cc.z())) - radsq;


  CubitBoolean inside_csph;
  if (dist2 > csTol) 
  {
    inside_csph = CUBIT_FALSE;
  }
  else if (dist2 < -csTol) 
  {
    inside_csph = CUBIT_TRUE;
  }
  else 
  {
    inside_csph = CUBIT_TRUE; // (numerically on the sphere)
  }

  if (inside_csph) 
  {
    // add it to the list and go check its neighbors

    neighbor_tet_list.append(tet_ptr);

    int iface;
    TET *adj_tet = NULL;
    NODE *na, *nb, *nc, *nd;
    tet_ptr->tet_nodes( na, nb, nc, nd );
    for (iface=0; iface<4; iface++) 
    {
      switch(iface)
      {
        case 0: adj_tet = (TET *)tet_ptr->get_connected_tet( nb, nc, nd ); break;
        case 1: adj_tet = (TET *)tet_ptr->get_connected_tet( nd, nc, na ); break;
        case 2: adj_tet = (TET *)tet_ptr->get_connected_tet( nd, na, nb ); break;
        case 3: adj_tet = (TET *)tet_ptr->get_connected_tet( nc, nb, na ); break;
      }
      if (adj_tet != NULL)
      {
        if (tet_visited(adj_tet) != tetVisited)
        {
          point_in_circumsphere( adj_tet, xx, neighbor_tet_list );
        }
      }
    }
  }
  return inside_csph;
}


//-------------------------------------------------------------------------
// Function:    circumsphere
// Description: get the circumsphere info for the tet
//              define a tooldata to hold the info if needed
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
circumsphere( TET *tet_ptr, CubitVector &center, double &radius_squared )
{
  ToolData *td = tet_ptr->get_TD( TDDelaunay< TET, NODE >::is_delaunay );
	TDDelaunay< TET, NODE > *td_del = dynamic_cast<TDDelaunay< TET, NODE >*> (td);
	if(td_del == NULL) {
		td_del = new TDDelaunay<TET, NODE>();
		tet_ptr->add_TD( td_del );
	}

	return td_del->circumsphere( tet_ptr, center, radius_squared );
}

//-------------------------------------------------------------------------
// Function:    set_tet_visited
// Description: set the visited flag for the tet
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__V
set_tet_visited( TET *tet_ptr, int new_visit_flag )
{
  ToolData *td = tet_ptr->get_TD( TDDelaunay< TET, NODE >::is_delaunay );
	TDDelaunay< TET, NODE > *td_del = dynamic_cast<TDDelaunay< TET, NODE >*> (td);
	if(td_del == NULL) {
		td_del = new TDDelaunay<TET, NODE>();
		tet_ptr->add_TD( td_del );
	}

	td_del->visit_flag( new_visit_flag );
}

//-------------------------------------------------------------------------
// Function:    tet_visited
// Description: return the visited flag for the tet
// Author:      sjowen
// Date:        8/2/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__I
tet_visited( TET *tet_ptr )
{
  ToolData *td = tet_ptr->get_TD( TDDelaunay< TET, NODE >::is_delaunay );
  if (td == NULL)
    return CUBIT_INT_MIN;
	TDDelaunay< TET, NODE > *td_del = dynamic_cast<TDDelaunay< TET, NODE >*> (td);
	if(td_del == NULL) 
    return CUBIT_INT_MIN;
  return td_del->visit_flag();
}

//-------------------------------------------------------------------------
// Function:    tet_volume
// Description: return the signed volume contained by the four coordinates
/*
                   top view

                       B
                       *
                      /|\
                     / | \
                    /  *D \
                   / _/ \_ \
                  /_/     \_\
               A *-----------* C
*/ 
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
TET_FACETOR_TOOL__D
tet_volume( const CubitVector &a, const CubitVector &b, 
            const CubitVector &c, const CubitVector &d )
{
   CubitVector da = a - d;
   CubitVector db = b - d;
   CubitVector dc = c - d;
   double vol = da.x()*(db.y()*dc.z() - dc.y()*db.z()) + 
                da.y()*(db.z()*dc.x() - dc.z()*db.x()) +
                da.z()*(db.x()*dc.y() - dc.x()*db.y());
   vol /= 6.0;

   return vol;
}

//-----------------------------------------------------------------------------
// Function: read_data
// Description: reads a "data" file -- a list of nodes and scalar/vector data
// Author: sjowen
// Date:8/13/2003
//------------------------------------------------------------------------------
TET_FACETOR_TOOL__I
read_data( const char *filename, DLIList<NODE *>&node_list )
{
  // open the file

  FILE *fp = fopen(filename, "r");
  if (fp == NULL)
  {
    PRINT_ERROR("Could not open file %s for reading\n", filename);
    return CUBIT_FAILURE;
  }

  PRINT_INFO("Reading data...\n");

  // read the number of nodes

  int iline = 1;
  char line[128];
  if (fgets(line, 128, fp) == NULL)
  {
    PRINT_ERROR("Format error in facet file %s on line %d\n", filename, iline);
    return CUBIT_FAILURE;
  }
  
  int nnodes = 0;
  int n = sscanf(line, "%d", &nnodes);
  if (n != 1)
  {
    PRINT_ERROR("Format error in data file %s on line %d\n", filename, iline);
    return CUBIT_FAILURE;
  }
  if (nnodes <= 0)
  {
    PRINT_ERROR("Expecting number of nodes in data file %s on line %d\n", filename, iline);
    return CUBIT_FAILURE;
  }

  // read the nodes

  int ii;
  double xx, yy, zz, data;
  NODE *new_node;
  for (ii=0; ii<nnodes; ii++)
  {
    iline++;
    if (fgets(line, 128, fp)== NULL)
    {
      PRINT_ERROR("Format error in data file %s on line %d\n", filename, iline);
      return CUBIT_FAILURE;
    }

    n = sscanf(line, "%lf %lf %lf %lf", &xx, &yy, &zz, &data );
    if (n < 3 || n > 4)
    {
      PRINT_ERROR("Format error in data file %s on line %d\n", filename, iline);
      PRINT_INFO("  Expecting 3 doubles and an optional data value, but instead read %d values\n", n);
      return CUBIT_FAILURE;
    }  
    CubitVector pt( xx, yy, zz );
    new_node = (NODE *) new NODE( pt );
    node_list.append( new_node );

    if (n==4)
    {
      TDInterpNode *td_interp = new TDInterpNode(data);
		  new_node->add_TD( td_interp );
    }
  }
  return CUBIT_SUCCESS;

}


// EOF

