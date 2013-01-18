//=======================================================================
//
//  File: CubitFacetEdgeData
//  Description: used for edges of CubitFacets.  These are optional
//               and used only if specific information must be stored at 
//               the edge common to tow facets
//  Owner: sjowen
//
//=======================================================================

#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitFacetData.hpp"
#include "CubitPoint.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "ToolData.hpp"

static int counter_id = 0;

//======================================================================
// Function: CubitFacetEdgeData (PUBLIC)
// Description: constructor.  Determines the adjacencies from existing
//              CubitPoint adjacency information
//              Note: can handle non-manifold topology
// Author: sjowen
// Date: 8/00
//======================================================================
CubitFacetEdgeData::CubitFacetEdgeData( CubitPoint *p1, CubitPoint *p2 )
{
  assert( p1 && p2 );
  assert( p1 != p2 );
  pointArray[0] = p1;
  pointArray[1] = p2;

  counter_id++;
  entityId = counter_id;
  
  // determine adjacency (assumes facets are already defined)
  // set the edge and use on the adjacent facets

  int sense;
  int eindex;
  CubitFacet *facet_ptr;
  DLIList<CubitFacet*> facet_list;
  p1->facets( facet_list );
  CubitBoolean found;
  for(int ii=0; ii<facet_list.size(); ii++) {
    found = CUBIT_FALSE;
    facet_ptr = facet_list.get_and_step();
    for(int jj=0; jj<3 && !found; jj++) 
    {
      if (facet_ptr->point(jj) == p2) 
      {
        found = CUBIT_TRUE;
        adjFacetList.append(facet_ptr);
        eindex = facet_ptr->edge_index( p1, p2, sense );
        facet_ptr->edge( this, eindex );
        facet_ptr->edge_use( sense, eindex );
      }
    }
  }
}

//======================================================================
// Function: CubitFacetEdgeData (PUBLIC)
// Description: overloaded constructor:  Adjacency information is
//              provided in the arguments -- assumes a 2D topology
// Author: sjowen
// Date: 8/00
//======================================================================
CubitFacetEdgeData::CubitFacetEdgeData( CubitPoint *p1, CubitPoint *p2,
                                CubitFacet *facet1, CubitFacet *facet2,
                                int edge_index1, int edge_index2)
{
  assert( p1 && p2 );
  assert( p1 != p2 );
  pointArray[0] = p1;
  pointArray[1] = p2;

  counter_id++;
  entityId = counter_id;
    
  if (facet1) {
    adjFacetList.append(facet1);
    facet1->edge( this, edge_index1 );
    facet1->edge_use( 1, edge_index1 );
  }
  if (facet2) {
    adjFacetList.append(facet2);
    facet2->edge( this, edge_index2 );
    facet2->edge_use( -1, edge_index2 );
  }
}

//======================================================================
// Function: ~CubitFacetEdgeData (PUBLIC)
// Description: destructor
// Author: sjowen
// Date: 8/00
//======================================================================
CubitFacetEdgeData::~CubitFacetEdgeData()
{
  assert(adjFacetList.size() == 0);
}

//======================================================================
// Function: adj_facet (PUBLIC)
// Description: get a specific adjacent facet (by index)
// Author: sjowen
// Date: 4/01
//======================================================================
CubitFacet *CubitFacetEdgeData::adj_facet( int index )
{
  if (index < 0 || index > adjFacetList.size() - 1)
    return (CubitFacet *)NULL;
  adjFacetList.reset();
  return adjFacetList.next(index);
}

//======================================================================
// Function: facets (PUBLIC)
// Description: get the atached facets
// Author: sjowen
// Date: 4/01
//======================================================================
void CubitFacetEdgeData::facets(DLIList<CubitFacet*> &facet_list)
{
  facet_list += adjFacetList;
}

//======================================================================
// Function: edges (PUBLIC)
// Description: get the atached edges
// Author: sjowen
// Date: 4/01
//======================================================================
void CubitFacetEdgeData::edges(DLIList<CubitFacetEdge*> &edge_list )
{
  edge_list.append( this );
}

//======================================================================
// Function: points (PUBLIC)
// Description: get the attached points
// Author: sjowen
// Date: 4/01
//======================================================================
void CubitFacetEdgeData::points(DLIList<CubitPoint*> &point_list )
{
  point_list.append( pointArray[0] );
  point_list.append( pointArray[1] );
}

//======================================================================
// Function: remove_facet (PUBLIC)
// Description: remove the facet from the edge adjacency list
// Author: sjowen
// Date: 4/01
//======================================================================
CubitStatus CubitFacetEdgeData::remove_facet( CubitFacet *facet_ptr )
{
  CubitStatus stat = CUBIT_SUCCESS;
  bool removed = adjFacetList.remove( facet_ptr );
  if (removed)
    stat = CUBIT_SUCCESS;
  else
    stat = CUBIT_FAILURE;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : Merge edges
//
// Special Notes : points must already be merged.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/24/03
//-------------------------------------------------------------------------
CubitStatus CubitFacetEdgeData::merge_edges( CubitFacetEdgeData* other_edge )
{
  if( other_edge == this )
    return CUBIT_SUCCESS;
  
  int reversed = false;
  if ( point(0) == other_edge->point(1) && point(1) == other_edge->point(0) )
    reversed = true;
  else if( point(0) != other_edge->point(0) || point(1) != other_edge->point(1) )
    return CUBIT_FAILURE;
  
  CubitFacetData* facet;
  adjFacetList.reset();
  while( other_edge->adjFacetList.size() )
  {
    facet = dynamic_cast<CubitFacetData*>(other_edge->adjFacetList.pop());
    for ( int i = 0; i < 3; i++ )
    {
      if( facet->edge(i) == other_edge )
      {
        facet->edge( this, i );
        add_facet(facet);
        if( reversed )
          facet->edge_use( -facet->edge_use(i), i );
      }
    }
  }
  
  delete other_edge;
  return CUBIT_SUCCESS;
}

void CubitFacetEdgeData::flip()
{
  CubitPoint* tmp = pointArray[0];
  pointArray[0] = pointArray[1];
  pointArray[1] = tmp;
  for ( int i = adjFacetList.size(); i--; )
  {
    CubitFacet* facet = adjFacetList.get_and_step();
    int index = facet->edge_index(this);
    assert(index >= 0);
    facet->edge_use( -facet->edge_use(index), index );
  }
  toggle_is_flipped();
}

