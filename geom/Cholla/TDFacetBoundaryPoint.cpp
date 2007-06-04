//- Class:       TDFacetBoundaryPoint
//- Description: Tool data for storing additional information at 
//-              the boundary of a facet set
//- Owner:       Steve Owen
//- Checked by:
//- Version:

#include "TDFacetBoundaryPoint.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitQuadFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CastTo.hpp"
#include "CubitTransformMatrix.hpp"
#include "GfxDebug.hpp"
TDFacetBoundaryPoint::TDFacetBoundaryPoint()
{

}

TDFacetBoundaryPoint::~TDFacetBoundaryPoint()
{
  for (int ii=0; ii<pointDataList.size(); ii++)
  {
    BoundaryPointData *bpd_ptr = pointDataList.get_and_step();
    delete bpd_ptr;
  }
}


//-------------------------------------------------------------------------
// Purpose       : create a new facet boundary edge
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::add_facet_boundary_point( 
  CubitPoint *point_ptr )
{
  ToolData *td;
  td = point_ptr->get_TD(&TDFacetBoundaryPoint::is_facet_boundary_point);
  if ( td == NULL )
  {
    TDFacetBoundaryPoint *td_gm = new TDFacetBoundaryPoint;
    point_ptr->add_TD( td_gm);
    td_gm->set_point( point_ptr );
  }
  else
  {
    TDFacetBoundaryPoint *td_gm = CAST_TO(td, TDFacetBoundaryPoint);
    td_gm->set_point( point_ptr );
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : create a new facet boundary point
//
// Special Notes : this is used to create and initialize with one facet
//                 and its corresponding normal.
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::add_facet_boundary_point( 
  CubitPoint *point_ptr,
  CubitFacet *facet_ptr,
  CubitVector &pt_normal )
{
  ToolData *td;
  TDFacetBoundaryPoint *td_gm = NULL;
  td = point_ptr->get_TD(&TDFacetBoundaryPoint::is_facet_boundary_point);
  if ( td == NULL )
  {
    td_gm = new TDFacetBoundaryPoint;
    point_ptr->add_TD( td_gm);
    td_gm->set_point( point_ptr );
  }
  else
  {
    td_gm = CAST_TO(td, TDFacetBoundaryPoint);
    td_gm->set_point( point_ptr );
  }

  BoundaryPointData *bpd_ptr = new BoundaryPointData;
  bpd_ptr->surfFacetList.append( facet_ptr ); 
  bpd_ptr->surfID = -1;
  bpd_ptr->uVal = 0.0; 
  bpd_ptr->vVal = 0.0;
  bpd_ptr->sizeVal = 0.0;
  bpd_ptr->normal = pt_normal;
  td_gm->pointDataList.append( bpd_ptr );
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : create a new facet boundary point
//
// Special Notes : this is used to create and initialize with one facet
//                 and its corresponding normal.  Same as above but uses 
//                 a CubitQuadFacet
//
// Creator       : Steve Owen
//
// Creation Date : 01/2004
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::add_facet_boundary_point( 
  CubitPoint *point_ptr,
  CubitQuadFacet *qfacet_ptr,
  CubitVector &pt_normal )
{
  ToolData *td;
  TDFacetBoundaryPoint *td_gm = NULL;
  td = point_ptr->get_TD(&TDFacetBoundaryPoint::is_facet_boundary_point);
  if ( td == NULL )
  {
    td_gm = new TDFacetBoundaryPoint;
    point_ptr->add_TD( td_gm);
    td_gm->set_point( point_ptr );
  }
  else
  {
    td_gm = CAST_TO(td, TDFacetBoundaryPoint);
    td_gm->set_point( point_ptr );
  }

  CubitFacet *facet0 = qfacet_ptr->get_tri_facet(0);
  CubitFacet *facet1 = qfacet_ptr->get_tri_facet(1);

  BoundaryPointData *bpd_ptr = new BoundaryPointData;
  if (facet0->point_index( point_ptr ) >= 0)
    bpd_ptr->surfFacetList.append( facet0 );
  if (facet1->point_index( point_ptr ) >= 0)
    bpd_ptr->surfFacetList.append( facet1 );
  bpd_ptr->surfID = -1;
  bpd_ptr->uVal = 0.0; 
  bpd_ptr->vVal = 0.0;
  bpd_ptr->sizeVal = 0.0;
  bpd_ptr->normal = pt_normal;
  td_gm->pointDataList.append( bpd_ptr );
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get the facet boundary point from a point
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
TDFacetBoundaryPoint* TDFacetBoundaryPoint::get_facet_boundary_point( 
  CubitPoint *point_ptr )
{
  ToolData *td;
  td = point_ptr->get_TD(&TDFacetBoundaryPoint::is_facet_boundary_point);
  if ( td != NULL )
  {
    TDFacetBoundaryPoint *td_gm = CAST_TO(td, TDFacetBoundaryPoint);
    return td_gm;
  }
  return (TDFacetBoundaryPoint*) NULL;
}


//-------------------------------------------------------------------------
// Purpose       : add a group of facets that are adjacent to this point
//                 that are on the same surface
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
void TDFacetBoundaryPoint::add_surf_facets(
  DLIList<CubitFacet *> adj_facet_list )
{
  BoundaryPointData *bpd_ptr = new BoundaryPointData;
  int ii;
  for (ii=0; ii<adj_facet_list.size(); ii++)
  {
    CubitFacet *facet_ptr = adj_facet_list.get_and_step();
    bpd_ptr->surfFacetList.append( facet_ptr );
  }
  
  bpd_ptr->surfID = -1;
  bpd_ptr->uVal = 0.0; 
  bpd_ptr->vVal = 0.0;
  bpd_ptr->sizeVal = 0.0;
  pointDataList.append( bpd_ptr );
  init_normal( bpd_ptr );
}

//-------------------------------------------------------------------------
// Purpose       : get the normal from a boundary point
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::get_normal( int surf_id, 
                                              CubitVector &normal )
{
  int found = 0;
  BoundaryPointData *bpd_ptr = NULL;
  for (int ii=0; ii<pointDataList.size(); ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    if(bpd_ptr->surfID == surf_id)
      found = 1;
  }
  if (!found)  // the surf_id does not match any in the list
    return CUBIT_FAILURE;

  normal = bpd_ptr->normal;

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get the normal from a boundary point based on its
//                 adjacent facet
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::get_normal( CubitFacet *adj_facet, 
                                              CubitVector &normal )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  normal = bpd_ptr->normal;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get the normal from a boundary point based on its
//                 adjacent edge
//
// Special Notes :  this is really only unique if the edge is non-manifold
//                  or the facets adjacent the edge all lie on the same 
//                  surface
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::get_normal( CubitFacetEdge *edge_ptr, 
                                              CubitVector &normal )
{
  DLIList<CubitFacet *>adj_facets;
  edge_ptr->facets(adj_facets);
  BoundaryPointData *bpd_ptr = NULL;
  for (int ii=0; ii<adj_facets.size() && bpd_ptr == NULL; ii++)
  {
    bpd_ptr = get_bpd( adj_facets.get_and_step() );
  }
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  normal = bpd_ptr->normal;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : set the normal into a boundary point
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::set_normal( int surf_id, 
                                              CubitVector &normal )
{
  int found = 0;
  BoundaryPointData *bpd_ptr = NULL;
  for (int ii=0; ii<pointDataList.size(); ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    if(bpd_ptr->surfID == surf_id)
      found = 1;
  }
  if (!found)  // the surf_id does not match any in the list
    return CUBIT_FAILURE;

  bpd_ptr->normal = normal;

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : reset the all normals for boundary points
//
// Special Notes :
//
// Creator       : Michael Brewer
//
// Creation Date : 02/05
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryPoint::reset_normals(  )
{
  BoundaryPointData *bpd_ptr;
  for (int ii=0; ii<pointDataList.size(); ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    if(bpd_ptr)
      init_normal(bpd_ptr);
    else{
      PRINT_ERROR("Could not determine boundary point data.\n");
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : set the surface id of associated with one of the 
//                 facets adjacent to this point
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
void TDFacetBoundaryPoint::set_surf_id( CubitFacet *facet_ptr, int surf_id )
{
  int found = 0;
  int ii, jj;
  CubitFacet *check_facet_ptr;
  BoundaryPointData *bpd_ptr;
  for (ii=0; ii<pointDataList.size() && !found; ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    for (jj=0; jj<bpd_ptr->surfFacetList.size() && !found; jj++)
    {
      check_facet_ptr = bpd_ptr->surfFacetList.get_and_step();
      if (check_facet_ptr == facet_ptr)
      {
        found = 1;
        bpd_ptr->surfID = surf_id;
      }
    }
  }
  assert(found);  // couldn't find the facet adjacent the surface
}

//-------------------------------------------------------------------------
// Purpose       : compute and set the normal
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
void TDFacetBoundaryPoint::init_normal( BoundaryPointData *bpd_ptr )
{
  
  if (bpd_ptr && bpd_ptr->surfFacetList.size() > 0) {
    CubitFacet *facet;
    double angle;
    int j;
    CubitVector normal;
    CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
    double totangle = 0.0e0;

      // weight the normal by the spanning angle at the point
    int mydebug =0;
    for (j = 0; j < bpd_ptr->surfFacetList.size(); j++)
    {
      facet = bpd_ptr->surfFacetList.get_and_step();
      angle = facet->angle( pointPtr );
      facet->weight( angle );
      totangle += angle;
      if(mydebug){
        if(angle <= 0.0){
          PRINT_INFO("Small angle.\n");
        }
      }
    }
      //First...
      // Attempt to handle this very odd case.  There is a facet with
      // zero area on the boundary.  Set normal to that of an adjacet
      // facet that doesn't have zero area.
    if(totangle == 0 && bpd_ptr->surfFacetList.size() == 1){
      PRINT_WARNING("Degenerate facet on edge of boundary.\n");
      facet = bpd_ptr->surfFacetList.get_and_step();
      
      CubitPoint* my_pt[3];
      my_pt[0] = facet->point(0);
      my_pt[1] = facet->point(1);
      my_pt[2] = facet->point(2);
      double largest_area=-1.0;
      double current_area=0.0;
      if(my_pt[0] == NULL || my_pt[1] == NULL || my_pt[2] == NULL){
        PRINT_ERROR("Problem determining normal.\n");
        return;
      }
      
      int tmp_i;
      for(tmp_i = 0; tmp_i < 3; tmp_i++){
        CubitFacet* other_facet = facet->shared_facet(my_pt[tmp_i],
                                                      my_pt[(tmp_i+1)%3]);
        if(other_facet != NULL){
          if(mydebug)
            other_facet->debug_draw(CUBIT_WHITE);
          current_area=other_facet->area();
          if(current_area>largest_area){
            largest_area=current_area;
            normal=other_facet->normal();
          }
        }
      }
      if(largest_area < CUBIT_DBL_MIN){
        PRINT_ERROR("Could not initialize facet normal.\n");
        return;
      }
      if(mydebug){
        facet->debug_draw(CUBIT_RED);
        PRINT_INFO("\n\nLargest area = %f\n\n",largest_area);
        GfxDebug::mouse_xforms();
      }
      bpd_ptr->normal = normal;
              
    }
      //Now...
      // Handle the normal case.
    else{
      for (j = 0; j < bpd_ptr->surfFacetList.size(); j++)
      {
        facet = bpd_ptr->surfFacetList.get_and_step();
        normal = facet->normal();
        normal.normalize();
        avg_normal += (facet->weight() / totangle) * normal;
      }
      avg_normal.normalize();
      bpd_ptr->normal = avg_normal;
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : return the boundary point data associated with a facet
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
BoundaryPointData *TDFacetBoundaryPoint::get_bpd( CubitFacet *facet )
{
  BoundaryPointData *bpd_ptr = NULL;
  int found = 0;
  int ii, jj;
  for (ii=0; ii<pointDataList.size() && !found; ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    for (jj=0; jj<bpd_ptr->surfFacetList.size() && !found; jj++)
    {
      if (bpd_ptr->surfFacetList.get_and_step() == facet)
        return bpd_ptr;
    }
  }
  
  return NULL;
}

//===============================================================================
//Function:  merge_normals (PRIVATE)
//Description: Enforce continuity across facets by merging the normals at two
//             points.  Assumes the points are coincident.
//Date: 05/01
//===============================================================================
CubitStatus TDFacetBoundaryPoint::merge_normals( 
  CubitFacet *facet0,
  CubitFacet *facet1)
{

  BoundaryPointData *bpd0_ptr = get_bpd( facet0 );
  BoundaryPointData *bpd1_ptr = get_bpd( facet1 );
  if (!bpd0_ptr || !bpd1_ptr)
    return CUBIT_FAILURE;

  CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
  double totangle = 0.0e0;
  double angle;
  CubitFacet *facet;
  int j;

  // weight the normal by the spanning angle at the point

  for (j = 0; j < bpd0_ptr->surfFacetList.size(); j++)
  {
    facet = bpd0_ptr->surfFacetList.get_and_step();
    angle = facet->angle( pointPtr );
    facet->weight( angle );
    totangle += angle;
  }

  for (j = 0; j < bpd1_ptr->surfFacetList.size(); j++)
  {
    facet = bpd1_ptr->surfFacetList.get_and_step();
    angle = facet->angle( pointPtr );
    facet->weight( angle );
    totangle += angle;
  }

  // computed weighted normal

  CubitVector normal;
  for (j = 0; j < bpd0_ptr->surfFacetList.size(); j++)
  {
    facet = bpd0_ptr->surfFacetList.get_and_step();
    normal = facet->normal();
    normal.normalize();
    avg_normal += (facet->weight() / totangle) * normal;
  }
  
  // orientation of facets may be opposite on opposing surfaces.
  // Check for this case and correct of necessary
  
  CubitVector norm0, norm1;
  norm0 = bpd0_ptr->normal;
  norm0.normalize();
  norm1 = bpd1_ptr->normal;
  norm1.normalize();
  double dot = norm0 % norm1;
  double sign = 1.0;
  if (dot < 0.0)
    sign = -1.0;

  for (j = 0; j < bpd1_ptr->surfFacetList.size(); j++)
  {
    facet = bpd1_ptr->surfFacetList.get_and_step();
    normal = sign * facet->normal();
    normal.normalize();
    avg_normal += (facet->weight() / totangle) * normal;
  }

  // set the new normal on both points

  avg_normal.normalize();
  bpd0_ptr->normal = avg_normal;
  CubitVector temp_vector = sign * avg_normal;
  bpd1_ptr->normal = temp_vector;

  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  set_uv (PUBLIC)
//Description: set the u-v values on a boundary point data
//Date: 06/01
//===============================================================================
CubitStatus TDFacetBoundaryPoint::set_uv( CubitFacet *adj_facet, 
                                          double u, double v )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet ); 
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  bpd_ptr->uVal = u;
  bpd_ptr->vVal = v;
  return CUBIT_SUCCESS;
}
//===============================================================================
//Function:  set_uvs (PUBLIC)
//Description: set the u-v values and size on a boundary point data
//Date: 06/01
//===============================================================================
CubitStatus TDFacetBoundaryPoint::set_uvs( CubitFacet *adj_facet, 
                                          double u, double v, double s )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet ); 
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  bpd_ptr->uVal = u;
  bpd_ptr->vVal = v;
  bpd_ptr->sizeVal = s;
  return CUBIT_SUCCESS;
}
//===============================================================================
//Function:  u (PUBLIC)
//Description: get the u value from a boundary point data
//Date: 06/01
//===============================================================================
double TDFacetBoundaryPoint::u( CubitFacet *adj_facet )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );                   
  assert( bpd_ptr != 0 );  // adjacent facet was not found in list
  return bpd_ptr->uVal;
}

//===============================================================================
//Function:  v (PUBLIC)
//Description: get the v value from a boundary point data
//Date: 06/01
//===============================================================================
double TDFacetBoundaryPoint::v( CubitFacet *adj_facet )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );                   
  assert( bpd_ptr != 0 );  // adjacent facet was not found in list
  return bpd_ptr->vVal;
}
//===============================================================================
//Function:  s (PUBLIC)
//Description: get the s value from a boundary point data
//Date: 07/02
//===============================================================================
double TDFacetBoundaryPoint::s( CubitFacet *adj_facet )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );                   
  assert( bpd_ptr != 0 );  // adjacent facet was not found in list
  return bpd_ptr->sizeVal;
}
//===============================================================================
//Function:  get_uv (PUBLIC)
//Description: get the uv values from a boundary point data
//Date: 06/01
//===============================================================================
CubitStatus TDFacetBoundaryPoint::get_uv( CubitFacet *adj_facet, 
                                          double &u, double &v )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );                   
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  u = bpd_ptr->uVal;
  v = bpd_ptr->vVal;
  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  get_uvs (PUBLIC)
//Description: get the uv values and size from a boundary point data
//Date: 7/02
//===============================================================================
CubitStatus TDFacetBoundaryPoint::get_uvs( CubitFacet *adj_facet, 
                                          double &u, double &v, double &s )
{
  BoundaryPointData *bpd_ptr = get_bpd( adj_facet );                   
  if (!bpd_ptr)
    return CUBIT_FAILURE;
  u = bpd_ptr->uVal;
  v = bpd_ptr->vVal;
  s = bpd_ptr->sizeVal;
  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  rotate_normal (PUBLIC)
//Description: rotate the normal at this point
//Date: 03/02
//===============================================================================
CubitStatus TDFacetBoundaryPoint::rotate_normal( CubitTransformMatrix &rotmat )
{
  BoundaryPointData * bpd_ptr;
  
  int ii;
  for (ii=0; ii<pointDataList.size(); ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    bpd_ptr->normal = rotmat * bpd_ptr->normal;
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  get_boundary_point_data_size (PUBLIC)
//Description: get the size of the boundary point data for mem alloc. purposes
//Date: 01/23/2003
//===============================================================================
CubitStatus TDFacetBoundaryPoint::get_boundary_point_data_size( 
   int &size_int_data,
   int &size_double_data )
{
   int num_bpd = pointDataList.size();
   int ii;
   BoundaryPointData *bpd_ptr;

   //pointPtr->id() and pointDataList.size()
   size_int_data += 2;
   for (ii=0; ii<num_bpd; ii++)
   {
     size_int_data += 2;
     bpd_ptr = pointDataList.get_and_step();
     //num_surfaces
     size_int_data += bpd_ptr->surfFacetList.size();
   }
   size_double_data += num_bpd * 6;
   return CUBIT_SUCCESS;
}
                                                   
//===============================================================================
//Function:  get_boundary_point_data (PUBLIC)
//Description: retreive the boundary point data so it can be dumped to a file
//             data must be allocated prior to calling this function!
//Date: 01/23/2003
//===============================================================================
CubitStatus TDFacetBoundaryPoint::get_boundary_point_data(int *int_data,
                                                          double *double_data,
                                                          int &iidx,
                                                          int &didx )
{
  int_data[iidx++] = pointPtr->id();
  int num_bpd = pointDataList.size();
  int_data[iidx++] = num_bpd;
  
  BoundaryPointData *bpd_ptr;
  int ii, jj;
  int numfacs = 0;
  CubitFacet *facet_ptr;
  for (ii=0; ii<num_bpd; ii++)
  {
    bpd_ptr = pointDataList.get_and_step();
    numfacs = bpd_ptr->surfFacetList.size();
    int_data[iidx++] = numfacs;
    for (jj=0; jj<bpd_ptr->surfFacetList.size(); jj++)
    {
      facet_ptr = bpd_ptr->surfFacetList.get_and_step();
      int_data[iidx++] = facet_ptr->id();
    }
    int_data[iidx++] = bpd_ptr->surfID;
    
    double_data[didx++] = bpd_ptr->normal.x();
    double_data[didx++] = bpd_ptr->normal.y();
    double_data[didx++] = bpd_ptr->normal.z();
    double_data[didx++] = bpd_ptr->sizeVal;
    double_data[didx++] = bpd_ptr->uVal;
    double_data[didx++] = bpd_ptr->vVal;
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  new_boundary_point_data (static PUBLIC)
//Description: create a new boundary point data and assign to a cubit point
//Author: sjowen
//Date: 01/26/2003
//===============================================================================
CubitStatus TDFacetBoundaryPoint::new_facet_boundary_point(CubitPoint **points,
                                                           CubitFacet **facets,
                                                           int &iidx,
                                                           int &didx,
                                                           int *int_data,
                                                           double *double_data)
{
  int id = int_data[iidx++];
  CubitPoint *point_ptr = points[id];
  TDFacetBoundaryPoint::add_facet_boundary_point(point_ptr);
  TDFacetBoundaryPoint *td = (TDFacetBoundaryPoint *)
    point_ptr->get_TD( &TDFacetBoundaryPoint::is_facet_boundary_point);

  td->initialize( facets, iidx, didx, int_data, double_data );

  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  initialize (PUBLIC)
//Description: initialize a boundary point data from data read from a CUB file
//Author: sjowen
//Date: 01/26/2003
//===============================================================================
void TDFacetBoundaryPoint::initialize(CubitFacet **facets,
                                      int &iidx,
                                      int &didx,
                                      int *int_data,
                                      double *double_data)
{
  int num_bpd = int_data[iidx++];
  
  BoundaryPointData *bpd_ptr;
  int ii, jj, id;
  int numfacs = 0;
  CubitFacet *facet_ptr;
  for (ii=0; ii<num_bpd; ii++)
  {
    bpd_ptr = new BoundaryPointData;
    pointDataList.append(bpd_ptr);
    numfacs = int_data[iidx++];
    for (jj=0; jj<numfacs; jj++)
    {
      id = int_data[iidx++];
      facet_ptr = facets[id];
      bpd_ptr->surfFacetList.append(facet_ptr);
    }
    bpd_ptr->surfID = int_data[iidx++];
    
    bpd_ptr->normal.x(double_data[didx++]);
    bpd_ptr->normal.y(double_data[didx++]);
    bpd_ptr->normal.z(double_data[didx++]);
    bpd_ptr->sizeVal = double_data[didx++];
    bpd_ptr->uVal = double_data[didx++];
    bpd_ptr->vVal = double_data[didx++];
  }

}

// eof
