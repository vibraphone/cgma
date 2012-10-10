#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "ToolData.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CastTo.hpp"

static int point_counter = 0;

//===========================================================================
// Function Name: CubitPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
CubitPointData::CubitPointData( double x_val, double y_val, double z_val )
  : coords( x_val, y_val, z_val )
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
}

//===========================================================================
// Function Name: CubitPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
CubitPointData::CubitPointData(double x_val, double y_val, double z_val,int *)
  : coords(x_val, y_val, z_val)
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
}

//===========================================================================
// Function Name: CubitPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
CubitPointData::CubitPointData( const CubitVector &new_point )
  : coords( new_point )
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
}

//===========================================================================
// Function Name: ~CubitPoint
//
// Member Type:  PUBLIC
// Description:  destructor
// Author: 
// Date:
//===========================================================================
CubitPointData::~CubitPointData()
{
  delete attachedFacets;
}

//===========================================================================
// Function Name: coordinates
//
// Member Type:  PUBLIC
// Description:  return coordinates of point in form of array
// Author: 
// Date:
//===========================================================================
void CubitPointData::coordinates(double point_array[3])
{
  point_array[0] = coords.x();
  point_array[1] = coords.y();
  point_array[2] = coords.z();
}

//===========================================================================
// Function Name: add_facet
//
// Member Type:  PUBLIC
// Description:  add a facet to the point's adjacency list
// Author: 
// Date:
//===========================================================================
void CubitPointData::add_facet( CubitFacet *facet )
{
  if ( attachedFacets == NULL )
    attachedFacets = new DLIList<CubitFacet*>(8);
  attachedFacets->append(facet);
  return;
}

//===========================================================================
// Function Name: remove_facet
//
// Member Type:  PUBLIC
// Description:  remove a facet to the point's adjacency list
// Author: 
// Date:
//===========================================================================
void CubitPointData::remove_facet( CubitFacet *facet )
{
  if ( attachedFacets == NULL )
    return;
  attachedFacets->remove(facet);
  return;
}

//===========================================================================
// Function Name: edges
//
// Member Type:  PUBLIC
// Description:  return the list of facet-edges attached to this point
// Author: sjowen
// Date:  4/28/01
//===========================================================================
void CubitPointData::edges( DLIList<CubitFacetEdge *> &edge_list )
{
  int ii, jj, kk;
  CubitFacet *facet_ptr;
  CubitFacetEdge *edge_ptr, *check_edge_ptr;
  
  if(attachedFacets == NULL)
    return;

  for (ii=0; ii<attachedFacets->size(); ii++)
  {
    facet_ptr = attachedFacets->get_and_step();
    for (jj=0; jj<3; jj++)
    {
      edge_ptr = facet_ptr->edge( jj );
      if (edge_ptr)
      {
        if (edge_ptr->point(0) == this ||
          edge_ptr->point(1) == this)
        {
          int found = 0;
          for (kk=0; kk<edge_list.size() && !found; kk++)
          {
            check_edge_ptr = edge_list.get_and_step();
            if (check_edge_ptr == edge_ptr)
              found = 1;
          }
          if (!found)
          {
            edge_list.append( edge_ptr );
          }
        }
      }
    }
  }
  
  return;
}


//===========================================================================
// Function Name: num_adj_facets
//
// Member Type:  PUBLIC
// Description:  return the number of facets attached to the point
// Author: 
// Date:
//===========================================================================
int CubitPointData::num_adj_facets()
{ 
  if (attachedFacets == NULL)
    return 0;
  else
    return attachedFacets->size(); 
}

//===========================================================================
// Function Name: merge_points
//
// Member Type:  PUBLIC
// Description:  merge this point with another
// Note:         currently ignores edges - call before defining the edges
// Author: sjowen
// Date: 9/18/01
//===========================================================================
CubitStatus CubitPointData::merge_points( CubitPoint *other_point,
                                          CubitBoolean keep_point /* = CUBIT_FALSE*/)
{
  if( other_point == this )
    return CUBIT_SUCCESS;
  
  DLIList<CubitFacet *>facet_list;
  other_point->facets( facet_list );
  CubitFacet *adj_facet;
  CubitFacetData *afd;
  int ii;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    adj_facet = facet_list.get_and_step();
    other_point->remove_facet( adj_facet );
    afd = CAST_TO(adj_facet, CubitFacetData);
    if (afd->point(0) == other_point)
      afd->set_point( this, 0 );
    else if(afd->point(1) == other_point)
      afd->set_point( this, 1 );
    else if(afd->point(2) == other_point)
      afd->set_point( this, 2 );
    else
    {
      assert(0);
      return CUBIT_FAILURE;
    }

    this->add_facet( adj_facet );
  
      // added by J.Kraftcheck - 2/14/03 - update edges also!!!
    CubitFacetEdgeData* afed;
    for ( int j = 0; j < 3; j++ ) {
      if( adj_facet->edge(j) ) {
        afed = dynamic_cast<CubitFacetEdgeData*>(adj_facet->edge(j));
        if ( afed->point(0) == other_point ) 
          afed->set_point(this, 0);
        else if( afed->point(1) == other_point )
          afed->set_point(this, 1);
      }
    } 
  }
  
  if (!keep_point)
    delete other_point;
  return CUBIT_SUCCESS;
}

CubitStatus CubitPointData::collapse_edge( CubitPointData *dead_point )
{
  int i, j;
  
    // Get the list of facets that will be destroyed
  DLIList<CubitFacet*> dead_facets, edge_facets;
  shared_facets( dead_point, dead_facets );

    // Get edges to update
  DLIList<CubitFacetEdge*> adj_edge_list;
  dead_point->edges(adj_edge_list);
  
    // Get edge to collapse
  CubitFacetEdge* collapse = 0;
  for ( i = adj_edge_list.size(); i-- && !collapse; )
    if (adj_edge_list.step_and_get()->other_point(dead_point) == this)
      collapse = adj_edge_list.extract();

    // For each dead facet...
  CubitFacetData *dead_facet, *edge_facet;
  CubitFacetEdgeData *keep_edge, *dead_edge;
  for ( i = dead_facets.size(); i--; )
  {
      // Get dead facet and relevant indices
    dead_facet = dynamic_cast<CubitFacetData*>(dead_facets.get_and_step());
    assert(!!dead_facet);
    int dead_pt_index = dead_facet->point_index(dead_point);
    int this_pt_index = dead_facet->point_index(this);
    int other_pt_index = 3 - dead_pt_index - this_pt_index;
    
      // The get the other point (facet should have 
      // this point, the dead point, and one other).
    CubitPoint* other_pt = dead_facet->point(other_pt_index);
    
      // Get the edges to merge
    dead_edge = dynamic_cast<CubitFacetEdgeData*>(dead_facet->edge(this_pt_index));
    keep_edge = dynamic_cast<CubitFacetEdgeData*>(dead_facet->edge(dead_pt_index));
    edge_facets.clean_out();

    //propagate all the tds from dead edge to keep edge
    DLIList<ToolData*> tds;
    dead_edge->get_all_TDs(&tds);
    for (int i=0; i<tds.size(); i++)
    {
      ToolData* new_td = tds.get_and_step()->propogate(keep_edge);
      if (new_td)
        keep_edge->add_TD(new_td);
    }
    
      // Get the list of adjacent facets to be updated
    dead_point->shared_facets(other_pt, edge_facets);
    
      // Determine relative sense of edges
    int dead_edge_pt = dead_edge->point(0) == other_pt ? 0 : 1;
    int rel_edge_sense = keep_edge->point(dead_edge_pt) == other_pt ? 1 : -1;
    
      // Update edge in each adjacent facet (merge keep_edge with dead_edge)
    for ( j = edge_facets.size(); j--; )
    {
      edge_facet = dynamic_cast<CubitFacetData*>(edge_facets.get_and_step());
      if ( edge_facet == dead_facet )
        continue;
      
      int dead_index = edge_facet->point_index(dead_point);
      int edge_index = (dead_index+1) % 3;
      if ( edge_facet->point(edge_index) == other_pt )
        edge_index = (dead_index+2) % 3;
      
      if ( dead_edge )
      {
        assert(edge_facet->edge(edge_index) == dead_edge);
        dead_edge->remove_facet(edge_facet);
        edge_facet->edge( 0, edge_index );
      }
      if ( keep_edge )
      {
        assert(edge_facet->edge(edge_index) == 0);
        keep_edge->add_facet(edge_facet);
        edge_facet->edge( keep_edge, edge_index );
        int use = rel_edge_sense * edge_facet->edge_use(edge_index);
        edge_facet->edge_use( use, edge_index );
      }
    }
    
      // delete dead entities
    
    delete dead_facet;
    
    if (dead_edge)
    {
      assert(dead_edge->num_adj_facets() == 0);
      
      adj_edge_list.move_to(dead_edge);
      assert(adj_edge_list.get() == dead_edge);
      adj_edge_list.extract();
      
      delete dead_edge;
    }
  }
  
  while (adj_edge_list.size())
  {
    keep_edge = dynamic_cast<CubitFacetEdgeData*>(adj_edge_list.pop());
    assert(!!keep_edge);
    if (keep_edge->point(0) == dead_point)
      keep_edge->set_point(this, 0);
    else if(keep_edge->point(1) == dead_point)
      keep_edge->set_point(this, 1);
    else
      assert(0);
  }
  
  DLIList<CubitFacet*> adj_facet_list(dead_point->num_adj_facets());
  dead_point->facets(adj_facet_list);
  while (adj_facet_list.size())
  {
    CubitFacetData* facet = dynamic_cast<CubitFacetData*>(adj_facet_list.pop());
    assert(!!facet);
    
    int index = facet->point_index(dead_point);
    assert((unsigned)index < (unsigned)3 && facet->point(index) == dead_point);
    
    dead_point->remove_facet(facet);
    this->add_facet(facet);
    facet->set_point( this, index );
  }
  
  
    // delete dead entities
  if ( collapse )
  {
    assert( collapse->num_adj_facets() == 0 );
    delete collapse;
  }

  assert( dead_point->num_adj_facets() == 0 );
  delete dead_point;
  
  return CUBIT_SUCCESS;
}
  

