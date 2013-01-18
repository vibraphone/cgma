#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "CubitPlane.hpp"
#include "CubitFacetEdge.hpp"
#include "FaceterFacetData.hpp"
#include "CastTo.hpp"
#include "ToolData.hpp"
#include "GfxDebug.hpp"

#ifndef CUBIT_MAX 
#define CUBIT_MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef CUBIT_MIN 
#define CUBIT_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define min3(a,b,c) CUBIT_MIN((CUBIT_MIN((a),(b))),(c))
#define max3(a,b,c) CUBIT_MAX((CUBIT_MAX((a),(b))),(c))
static int counter_id = 0;

//===========================================================================
//Function Name: FaceterFacetData
//
//Member Type:  PUBLIC
//Description:  constructor
//===========================================================================
FaceterFacetData::FaceterFacetData( CubitPoint *p1, CubitPoint *p2,
                                    CubitPoint *p3 )
{
  assert( p1 && p2 && p3 );
  assert( p1 != p2 && p1 != p3 && p2 != p3 );
  myOwner = NULL;
  pointArray[0] = p1;
  pointArray[1] = p2;
  pointArray[2] = p3;
  p1->add_facet(this);
  p2->add_facet(this);
  p3->add_facet(this);

  edgeArray[0] = NULL;
  edgeArray[1] = NULL;
  edgeArray[2] = NULL;
  edgeUse[0] = edgeUse[1] = edgeUse[2] = 0;
  patchCtrlPts = NULL;

  //CubitVector temp = normal();
  //cachedNormal = new CubitVector(temp);
  plane();
  counter_id++;
  entityId = counter_id;

  // define the bounding box

  CubitVector bbox_min, bbox_max;
  bbox_min.x(min3(p1->x(),p2->x(),p3->x()));
  bbox_min.y(min3(p1->y(),p2->y(),p3->y()));
  bbox_min.z(min3(p1->z(),p2->z(),p3->z()));
  bbox_max.x(max3(p1->x(),p2->x(),p3->x()));
  bbox_max.y(max3(p1->y(),p2->y(),p3->y()));
  bbox_max.z(max3(p1->z(),p2->z(),p3->z()));
  bBox.reset(bbox_min,bbox_max);

}

//===========================================================================
//Function Name: ~FaceterFacetData
//
//Member Type:  PUBLIC
//Description:  destructor
//===========================================================================
FaceterFacetData::~FaceterFacetData()
{
}


//===========================================================================
//Function Name: closest_point
//
//Member Type:  PUBLIC
//Description:  return the closest point on plane defined by the facet
//===========================================================================
CubitStatus FaceterFacetData::closest_point(const CubitVector &point, 
                                          CubitVector &closest_point )
{
  if( cachedPlane )
  {
    closest_point = cachedPlane->project( point );
    return CUBIT_SUCCESS;
  }
  
  CubitVector normal_vec = normal();
  CubitVector P0 = pointArray[0]->coordinates();
  CubitVector point_to_P0 = point - P0;
  CubitVector P0_to_projection;
  //Store the normal squared the the dot product of the normal and vector to the point.
  double normal_length_sq = normal_vec.length_squared();
  //If the normal is zero, get out!
  if ( normal_length_sq < CUBIT_RESABS )
     return CUBIT_FAILURE;
  double point_to_P0_dot_norm = point_to_P0%normal_vec;
  //Use P0_to_projection as a temporary to store half the equation.
  P0_to_projection = (point_to_P0_dot_norm/normal_length_sq)*normal_vec;
  P0_to_projection = point_to_P0 - P0_to_projection;
  //Now we have the vector from P0, to the closest point.  Add the vector 
  //to P0, and you have the location!
  closest_point = P0 + P0_to_projection;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Local modification functions.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/25/00
//-------------------------------------------------------------------------
CubitPoint* FaceterFacetData::split_edge( CubitPoint* edge1_pt, 
                                    CubitPoint* edge2_pt,
                                    const CubitVector& position )
{
  CubitPointData* new_pt = new CubitPointData(position);
  
    // split triangles
  
  DLIList<CubitFacet*> facets;
  edge1_pt->shared_facets( edge2_pt, facets );
  
  facets.reset();
  for ( int i = facets.size(); i--; ) {

    CubitFacet* facet = facets.get_and_step();
    FaceterFacetData* facet_d = dynamic_cast<FaceterFacetData*>(facet);
    assert(!!facet_d);
   
 
      // fix up existing facet
    
    int pt2_index = facet->point_index( edge2_pt );
    bool edge_reversed = ( edge1_pt == facet->point( (pt2_index+1) % 3 ) );
    int edge_index = (pt2_index + 2 - edge_reversed) % 3;
    
    facet->update_plane();
    edge2_pt->remove_facet( facet );
    facet_d->set_point( new_pt, pt2_index );
    new_pt->add_facet( facet );
    

      // make new facet
      
    CubitPoint* other_pt = facet->point( edge_index );
    FaceterFacetData* new_facet;
    if ( edge_reversed )
      new_facet = new FaceterFacetData( other_pt, edge2_pt, new_pt );
    else
      new_facet = new FaceterFacetData( other_pt, new_pt, edge2_pt );
  }
       
  return new_pt; 
}

//-------------------------------------------------------------------------
// Purpose       : insert a point into the facet
//
// Special Notes : create two new facets and return them
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitPoint* FaceterFacetData::insert_point( const CubitVector& position,
                                          CubitFacet*& new_tri1,
                                          CubitFacet*& new_tri2 )  
{
  if( DEBUG_FLAG(110 ) ) 
  {
    GfxDebug::draw_facet(this, CUBIT_WHITE);
    GfxDebug::flush();
  }
  CubitPoint* new_point = (CubitPoint *) new CubitPointData( position );
  
  new_tri1 = (CubitFacet *)
    new FaceterFacetData( pointArray[1], pointArray[2], new_point );
  new_tri2 = (CubitFacet *)
    new FaceterFacetData( pointArray[2], pointArray[0], new_point );
  
  pointArray[2]->remove_facet( this );
  pointArray[2] = new_point;
  new_point->add_facet( this );

  if( DEBUG_FLAG(110) )
  {
    CubitVector temp_vector = new_point->coordinates();
    GfxDebug::draw_facet(this, CUBIT_MAGENTA);
    GfxDebug::flush();
    GfxDebug::draw_facet(new_tri1, CUBIT_MAGENTA);
    GfxDebug::flush();
    GfxDebug::draw_facet(new_tri2, CUBIT_MAGENTA);
    GfxDebug::flush();
    GfxDebug::draw_point(temp_vector, CUBIT_CYAN);
    GfxDebug::flush();
  }
  
  if( cachedPlane && 
      (cachedPlane->distance(position) > CUBIT_RESABS) )
  {
    delete cachedPlane;
    cachedPlane = 0;
  }
  
  if (isBackwards)
  {
    new_tri1->is_backwards( 1 );
    new_tri2->is_backwards( 1 );
  }
  return new_point;
}


//-------------------------------------------------------------------------
// Purpose       : reorient the facet
//
// Special Notes :  This function should only be used if it is acceptable 
//                  to change the underlying data representation of the 
//                  facets  (for example the Sierra facet representation 
//                  cannot change and will not permit changing of the
//                  node orders on a facet) 
//
//                  Also note that if this function is used, the isBackwards
//                  flag on the facet should be set appropriately.
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void FaceterFacetData::flip()
{
  CubitPoint *pt_tmp = pointArray[1];
  pointArray[1] = pointArray[2];
  pointArray[2] = pt_tmp;

  CubitFacetEdge *ed_tmp = edgeArray[1];
  edgeArray[1] = edgeArray[2];
  edgeArray[2] = ed_tmp;

  for (int ii=0; ii<3; ii++)
  {
    if (edgeUse[ii] == -1)
      edgeUse[ii] = 1;
    else if(edgeUse[ii] == 1)
      edgeUse[ii] = -1;
  }
  if (cachedPlane)
  {
    CubitVector normal = cachedPlane->normal();
    normal = -normal;
    cachedPlane->normal( normal );
  }
}


