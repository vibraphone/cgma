#include "ImprintPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "ToolData.hpp"
#include "CubitFacetData.hpp"
#include "CastTo.hpp"
#include "RefEntity.hpp"
#include "CubitMessage.hpp"

static int point_counter = 0;

//===========================================================================
// Function Name: ImprintPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
ImprintPointData::ImprintPointData( double x_val, double y_val, double z_val )
  : coords( x_val, y_val, z_val )
{
  initialize_imprint_data();
}
void ImprintPointData::initialize_imprint_data()
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
  myOwner = NULL;
  myPointType = UNSET_POINT;
  ownerRefVertex = CUBIT_FALSE;
  matchingPoint = NULL;
  startPartition = CUBIT_FALSE;
  endPartition = CUBIT_FALSE;
  listLoopPos = 0;
  loopPos = 0;
  pointMatches = NULL;
  isMatched = CUBIT_FALSE;
  prevSeg = NULL;
  nextSeg = NULL;
}


  
CubitStatus ImprintPointData::set_matching_point(ImprintPointData *other)
{
 // assert(matchingPoint == NULL|| other == matchingPoint);
  matchingPoint = other;
  return CUBIT_SUCCESS;
}
//===========================================================================
// Function Name: ImprintPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
ImprintPointData::ImprintPointData( const CubitVector &new_point )
  : coords( new_point )
{
  initialize_imprint_data();
}

//===========================================================================
// Function Name: ImprintPointData
//
// Member Type:  PUBLIC
// Description:  copy constructor
// Author: 
// Date:
//===========================================================================
ImprintPointData::ImprintPointData( ImprintPointData *copy, CubitVector &new_coords )
{
  coords = new_coords;
  initialize_imprint_data();
  myPointType=copy->myPointType;
  owner(copy->myOwner);
  set_matching_point(copy);
  copy->set_matching_point(this);
}


//===========================================================================
// Function Name: ~CubitPoint
//
// Member Type:  PUBLIC
// Description:  destructor
// Author: 
// Date:
//===========================================================================
ImprintPointData::~ImprintPointData()
{
  if ( attachedFacets )
    delete attachedFacets;
  if ( pointMatches )
    delete pointMatches;
}

//===========================================================================
// Function Name: coordinates
//
// Member Type:  PUBLIC
// Description:  return coordinates of point in form of array
// Author: 
// Date:
//===========================================================================
void ImprintPointData::coordinates(double point_array[3])
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
void ImprintPointData::add_facet( CubitFacet *facet )
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
void ImprintPointData::remove_facet( CubitFacet *facet )
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
void ImprintPointData::edges( DLIList<CubitFacetEdge *> &edge_list)
{
  int ii, jj, kk;
  CubitFacet *facet_ptr;
  CubitFacetEdge *edge_ptr, *check_edge_ptr;
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
int ImprintPointData::num_adj_facets()
{ 
  if (attachedFacets == NULL)
    return 0;
  else
    return attachedFacets->size(); 
}

void ImprintPointData::add_match_data(ImprintMatchData *match_data)
{
  if ( !pointMatches )
    pointMatches = new DLIList<ImprintMatchData*>;
  pointMatches->append(match_data);
}
void ImprintPointData::get_match_list(DLIList <ImprintMatchData*> &match_data_list)
{
  if ( pointMatches != NULL )
  {
    match_data_list += *pointMatches;
  }
}
void ImprintPointData::set_unmatched()
{
  if ( pointMatches != NULL)
  {
    pointMatches->clean_out();
  }
  isMatched = CUBIT_FALSE;
}
void ImprintPointData::set_matched(ImprintMatchData *match_data)
{
  if ( !pointMatches )
    pointMatches = new DLIList<ImprintMatchData*>;
  else
  {
    pointMatches->clean_out();
  }
  pointMatches->append(match_data);
  isMatched = CUBIT_TRUE;
}
