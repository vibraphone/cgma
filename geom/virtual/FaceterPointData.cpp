#include "FaceterPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "ToolData.hpp"
#include "CubitFacetData.hpp"
#include "CastTo.hpp"
#include "RefEntity.hpp"

static int point_counter = 0;

//===========================================================================
// Function Name: FaceterPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
FaceterPointData::FaceterPointData( double x_val, double y_val, double z_val )
  : coords( x_val, y_val, z_val )
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
  myOwner = NULL;
  myInteriorAngle = CUBIT_DBL_MAX;
  mPrev = NULL;
  mNext = NULL;
}

//===========================================================================
// Function Name: FaceterPointData
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
FaceterPointData::FaceterPointData( const CubitVector &new_point )
  : coords( new_point )
{
  attachedFacets = NULL;
  point_counter++;
  entityId = point_counter;
  myOwner = NULL;
  myInteriorAngle = CUBIT_DBL_MAX;
  mPrev = NULL;
  mNext = NULL;
}

//===========================================================================
// Function Name: ~CubitPoint
//
// Member Type:  PUBLIC
// Description:  destructor
// Author: 
// Date:
//===========================================================================
FaceterPointData::~FaceterPointData()
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
void FaceterPointData::coordinates(double point_array[3])
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
void FaceterPointData::add_facet( CubitFacet *facet )
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
void FaceterPointData::remove_facet( CubitFacet *facet )
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
void FaceterPointData::edges( DLIList<CubitFacetEdge *> &edge_list)
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
int FaceterPointData::num_adj_facets()
{ 
  if (attachedFacets == NULL)
    return 0;
  else
    return attachedFacets->size(); 
}

int FaceterPointData::sort_by_angle (FaceterPointData *&pt_1,
                                     FaceterPointData *&pt_2)
{
  if ( (pt_1->get_interior_angle() > pt_2->get_interior_angle() ) )
    return -1;
  if ( (pt_2->get_interior_angle() < pt_2->get_interior_angle() ) )
    return 1;
  else
    return 0;
}


