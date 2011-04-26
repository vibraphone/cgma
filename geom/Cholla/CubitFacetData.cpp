#include "CubitFacet.hpp" 
#include "CubitPoint.hpp" 
#include "CubitPointData.hpp" 
#include "CubitVector.hpp" 
#include "GeometryDefines.h"  
#include "CubitPlane.hpp" 
#include "CubitFacetEdge.hpp" 
#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp" 
#include "CastTo.hpp" 
#include "ToolData.hpp" 
#include "GfxDebug.hpp"
#include "TDFacetBoundaryPoint.hpp"

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
//Function Name: CubitFacetData 
// 
//Member Type:  PUBLIC 
//Description:  constructor 
//=========================================================================== 
CubitFacetData::CubitFacetData( CubitPoint *p1, CubitPoint *p2, 
                                CubitPoint *p3 ) 
{ 
  assert( p1 && p2 && p3 ); 
  assert( p1 != p2 && p1 != p3 && p2 != p3 ); 
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
 
  plane(); 
  counter_id++; 
  entityId = counter_id; 
 
  define_bounding_box();
 
} 

//=========================================================================== 
//Function Name: CubitFacetData 
// 
//Member Type:  PUBLIC 
//Description:  constructor 
//=========================================================================== 
CubitFacetData::CubitFacetData( CubitFacetEdge *e1, CubitFacetEdge *e2, 
                                CubitFacetEdge *e3 ) 
{ 
  assert( e1 && e2 && e3 ); 
  assert( e1 != e2 && e1 != e3 && e2 != e3 ); 

  edgeArray[0] = e1; 
  edgeArray[1] = e2; 
  edgeArray[2] = e3; 

  define_point(e2, e3, 0);
  define_point(e3, e1, 1);
  define_point(e1, e2, 2);
  e1->add_facet(this);
  e2->add_facet(this);
  e3->add_facet(this);

  patchCtrlPts = NULL; 

  plane(); 
  counter_id++; 
  entityId = counter_id; 
 
  define_bounding_box();
 
} 

//=========================================================================== 
//Function Name: CubitFacetData 
// 
//Member Type:  PUBLIC 
//Description:  constructor 
//=========================================================================== 
CubitFacetData::CubitFacetData( CubitPoint *p1, CubitPoint *p2,
                                CubitPoint *p3, int *tool_data)
{
  assert( p1 && p2 && p3 );
  assert( p1 != p2 && p1 != p3 && p2 != p3 );
  pointArray[0] = p1;
  pointArray[1] = p2;
  pointArray[2] = p3;

  edgeArray[0] = NULL; 
  edgeArray[1] = NULL; 
  edgeArray[2] = NULL; 

  p1->add_facet(this);
  p2->add_facet(this);
  p3->add_facet(this);

  edgeArray[0] = 0;
  edgeArray[1] = 0;
  edgeArray[2] = 0;
  allocate_edge(p2,p3,0);
  allocate_edge(p3,p1,1);
  allocate_edge(p1,p2,2);

  patchCtrlPts = NULL;

  plane();
  counter_id++;
  entityId = counter_id;
  // update toolID
  set_tool_id(*tool_data);

  define_bounding_box();
} 
 
//=========================================================================== 
//Function Name: ~CubitFacetData 
// 
//Member Type:  PUBLIC 
//Description:  destructor 
//=========================================================================== 
CubitFacetData::~CubitFacetData() 
{
  destruct_facet_internals();
}

void CubitFacetData::destruct_facet_internals()
{
  int ii = 3; 
  for (ii = 2; ii>=0; ii--){
    //remove this triangle-point association at the points. 
    CubitPoint *current_point = point(ii);
    if (current_point)
      current_point->remove_facet(this);     

    pointArray[ii] = NULL;
    
     //remove edge-point association at the edges
     
     CubitFacetEdge *current_edge = edgeArray[ii];
     
     CubitStatus status;
     if (current_edge)
       status = current_edge->remove_facet(this);
     edgeArray[ii] = NULL;
  }
}


//=========================================================================== 
//Function Name: closest_point 
// 
//Member Type:  PUBLIC 
//Description:  return the closest point on plane defined by the facet 
//=========================================================================== 
CubitStatus CubitFacetData::closest_point(const CubitVector &point,  
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
CubitPoint* CubitFacetData::split_edge( int edge_index,
                                        const CubitVector& position )
{
  CubitPoint* pt1 = point((edge_index+1)%2);
  CubitPoint* pt2 = point((edge_index+2)%2);
  return split_edge( pt1, pt2, position );
}

CubitPoint* CubitFacetData::split_edge( CubitPoint* edge1_pt,  
                                        CubitPoint* edge2_pt, 
                                        const CubitVector& position ) 
{ 
  CubitPointData* new_pt = new CubitPointData(position);
  
    // split edge, if there is one
  
  CubitFacetEdge* edge = edge1_pt->shared_edge( edge2_pt );
  CubitFacetEdgeData* new_edge = 0;
  if ( edge ) {
    CubitFacetEdgeData* edge_d = dynamic_cast<CubitFacetEdgeData*>(edge);
    assert(!!edge_d);
    
      // make sure new edge has same orientation as old edge
    new_edge = dynamic_cast<CubitFacetEdgeData*>(new_pt->shared_edge(edge2_pt));
    if ( edge->point(0) == edge1_pt ) {
      edge_d->set_point(new_pt, 1);
      if ( !new_edge )
      {
        new_edge = new CubitFacetEdgeData( new_pt, edge2_pt );
        DLIList<ToolData*> tds;
        edge->get_all_TDs(&tds);
        for (int i=0; i<tds.size(); i++)
        {
          ToolData* new_td = tds.get_and_step()->propogate(new_edge);
          if (new_td)
            new_edge->add_TD(new_td);
        }
      }
      else if( new_edge->point(0) != new_pt )
        new_edge->flip();
    } else {
      edge_d->set_point(new_pt, 0);
      if ( !new_edge )
      {
        new_edge = new CubitFacetEdgeData( edge2_pt, new_pt );
        DLIList<ToolData*> tds;
        edge->get_all_TDs(&tds);
        for (int i=0; i<tds.size(); i++)
        {
          ToolData* new_td = tds.get_and_step()->propogate(new_edge);
          if (new_td)
            new_edge->add_TD(new_td);
        }
      }
      else if( new_edge->point(1) != new_pt )
        new_edge->flip();
    }
  }
  
    // split triangles
  
  DLIList<CubitFacet*> facets;
  edge1_pt->shared_facets( edge2_pt, facets );
  
  facets.reset();
  for ( int i = facets.size(); i--; ) {

    CubitFacet* facet = facets.get_and_step();
    CubitFacetData* facet_d = dynamic_cast<CubitFacetData*>(facet);
    assert(!!facet_d);
   
 
      // fix up existing facet
    
    int pt2_index = facet->point_index( edge2_pt );
    bool edge_reversed = ( edge1_pt == facet->point( (pt2_index+1) % 3 ) );
    int edge_index = (pt2_index + 1 + edge_reversed) % 3;
    
    edge2_pt->remove_facet( facet );
    facet_d->set_point( new_pt, pt2_index );
    new_pt->add_facet( facet );
    facet->update_plane();
    

      // make new facet
      
    CubitPoint* other_pt = facet->point( edge_index );
    CubitFacetData* new_facet;
    if ( edge_reversed )
      new_facet = new CubitFacetData( other_pt, edge2_pt, new_pt );
    else
      new_facet = new CubitFacetData( other_pt, new_pt, edge2_pt );

    DLIList<ToolData*> td_list;
    facet->get_all_TDs(&td_list);
    for (int i=0; i< td_list.size(); i++)
    {
      ToolData* new_td = td_list.get_and_step()->propogate(new_facet);
      if (new_td)
      {
        new_facet->add_TD(new_td);
      }
    }
   
    if ( new_edge ) {
      assert(!new_facet->edge(0));
      new_facet->edge( new_edge, 0 );
      new_edge->add_facet( new_facet );
      int sense = new_facet->point( 1 ) == new_edge->point(0) ? 1 : -1;
      new_facet->edge_use( sense, 0 );
    }

    
      // move other edge, if there is one

    int pt1_index =  (pt2_index + 2 - edge_reversed) % 3;
    CubitFacetEdge* other_edge = facet->edge(pt1_index);
    if ( other_edge ) {
      other_edge->remove_facet(facet);
      facet->edge( 0, pt1_index );
      int e_index = 1 + edge_reversed;
      assert(!new_facet->edge(e_index));
      new_facet->edge( other_edge, e_index );
      other_edge->add_facet(new_facet);
      int sense = new_facet->point( (e_index+1)%3 ) == other_edge->point(0) ? 1 : -1;
      new_facet->edge_use( sense, e_index );
    }

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
CubitPoint* CubitFacetData::insert_point( const CubitVector& position, 
                                          CubitFacet*& new_tri1_out, 
                                          CubitFacet*& new_tri2_out )   
{ 
  CubitPointData* new_point = new CubitPointData( position ); 
  CubitFacetData *new_tri1, *new_tri2;
  new_tri1 = new CubitFacetData( point(1), point(2), new_point ); 
  new_tri2 = new CubitFacetData( point(2), point(0), new_point ); 
   
  point(2)->remove_facet( this ); 
  set_point( new_point, 2 );
  new_point->add_facet( this );
  
  if ( edge(0) ) {
    new_tri1->edge( edge(0), 2 );
    new_tri1->edge_use( edge_use(0), 2 );
    edge(0)->remove_facet(this);
    edge(0)->add_facet(new_tri1);
    edge( 0, 0 );
  }
  
  if ( edge(1) ) {
    new_tri2->edge( edge(1), 2 );
    new_tri2->edge_use( edge_use(1), 2 );
    edge(1)->remove_facet(this);
    edge(1)->add_facet(new_tri2);
    edge( 0, 1 );
  }
 
  update_plane();
  
  new_tri1_out = new_tri1;
  new_tri2_out = new_tri2;
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
void CubitFacetData::flip() 
{ 
  CubitVector* ctrl_points=control_points( );
  
  CubitPoint *pt_tmp = pointArray[1]; 
  pointArray[1] = pointArray[2]; 
  pointArray[2] = pt_tmp; 
 
  CubitFacetEdge *ed_tmp = edgeArray[1]; 
  edgeArray[1] = edgeArray[2]; 
  edgeArray[2] = ed_tmp;

    //make sure the edgeUses are matched with the correct edge...
  int ed_use_tmp = edgeUse[1];
  edgeUse[1]=edgeUse[2];
  edgeUse[2]=ed_use_tmp;
  
    //now flip the edge uses...
  int ii;
  for (ii=0; ii<3; ii++) 
  { 
    if (edgeUse[ii] == -1) {
      edgeUse[ii] = 1; 
    }
    else if(edgeUse[ii] == 1) {
      edgeUse[ii] = -1; 
    }
  }
  if(ctrl_points){
    CubitVector tmp_point;
    tmp_point = ctrl_points[0];
    ctrl_points[0]=ctrl_points[1];
    ctrl_points[1]=tmp_point;
    tmp_point = ctrl_points[2];
    ctrl_points[2]=ctrl_points[5];
    ctrl_points[5]=tmp_point;
    tmp_point = ctrl_points[3];
    ctrl_points[3]=ctrl_points[4];
    ctrl_points[4]=tmp_point;
  }
  update_plane();
    //update the normals on the points (including boundary points)
  for (ii=0; ii<3; ii++) 
  {
    pointArray[ii]->compute_avg_normal();
    TDFacetBoundaryPoint* tdfbp =
      TDFacetBoundaryPoint::get_facet_boundary_point(pointArray[ii]);
    if(tdfbp){
      if(!tdfbp->reset_normals()){
        PRINT_ERROR("Could not reset all the normals for a point.\n");
      }
    }
  }
  
} 
 
//===========================================================================
//Function Name: allocate_edge
//
//Member Type:  PRIVATE
//Description:  associate facetedge with this facet
//Special Notes: 1. Assume edge(p1,p2) is in Counter-clockwise direction wrt to//                  this facet)
//               2. Helper function to Constructor
//===========================================================================
void CubitFacetData::allocate_edge(CubitPoint *p1, CubitPoint *p2, int edge_index){
  
  assert(edge_index >= 0 && edge_index < 3);
  assert(p1 != NULL && p2 != NULL);

  CubitFacetEdge* shared_edge = p1->get_edge(p2);
  
  if(shared_edge == NULL){
    //- if edge don't exist, create it
    edgeArray[edge_index] = (CubitFacetEdge *) new CubitFacetEdgeData(p1,p2);
    edgeUse[edge_index] = 1;
  } else {
    edgeArray[edge_index] = shared_edge;
    shared_edge->add_facet(this);
    if(shared_edge->point(0) == p1)
      edgeUse[edge_index]= 1;
    else
      edgeUse[edge_index] = -1;
    
  }
  
}


//===========================================================================
//Function Name: define_point
//
//Member Type:  PRIVATE
//Description:  define the point based on its adjacent edges.  define the 
//              edge uses too.
//Special Notes: 1. Assume edge(e1,e2) is in Counter-clockwise direction wrt 
//                  this facet
//               2. Helper function to Constructor
//===========================================================================
void CubitFacetData::define_point(CubitFacetEdge *e1, 
                                  CubitFacetEdge *e2, 
                                  int point_index)
{
  
  assert(point_index >= 0 && point_index < 3);
  assert(e1 != NULL && e2 != NULL);

  CubitPoint *pA, *pB, *pC, *pD;
  pA = e1->point(0);
  pB = e1->point(1);
  pC = e2->point(0);
  pD = e2->point(1);
  if(pC == pB || pC == pA)
  {
    pointArray[point_index] = pC;
    pC->add_facet(this);
    edgeUse[(point_index+2)%3] = 1;
  }
  else if(pD == pB || pD == pA)
  {
    pointArray[point_index] = pD;
    pD->add_facet(this);
    edgeUse[(point_index+2)%3] = -1;
  }
  else
  {
    assert(0);  // the edges are not adjacent;
  }
  
}

//===========================================================================
//Function Name: define_bounding_box
//
//Member Type:  PRIVATE
//Description:  compute and store the bounding box
//               Helper function to Constructor
//===========================================================================
void CubitFacetData::define_bounding_box()
{
  CubitVector bbox_min, bbox_max; 
  CubitPoint *p1 = pointArray[0];
  CubitPoint *p2 = pointArray[1];
  CubitPoint *p3 = pointArray[2];
  bbox_min.x(min3(p1->x(),p2->x(),p3->x())); 
  bbox_min.y(min3(p1->y(),p2->y(),p3->y())); 
  bbox_min.z(min3(p1->z(),p2->z(),p3->z())); 
  bbox_max.x(max3(p1->x(),p2->x(),p3->x())); 
  bbox_max.y(max3(p1->y(),p2->y(),p3->y())); 
  bbox_max.z(max3(p1->z(),p2->z(),p3->z())); 
  bBox.reset(bbox_min,bbox_max);
}

CubitStatus CubitFacetData::flip_edge( CubitFacetEdge *edge )
{
  int i;
  for(i=0; i<3; i++)
  {
    if (edgeArray[i] == edge)
      return flip_edge(i);
  }
  return CUBIT_FAILURE;
}

CubitStatus CubitFacetData::flip_edge( int this_edge_index ) 
{
    // get point indices on this facet
  int this_pt1_index = (this_edge_index+1)%3;
  int this_pt2_index = (this_edge_index+2)%3;
  
    // get edge points
  CubitPoint* point1 = point(this_pt1_index);
  CubitPoint* point2 = point(this_pt2_index);

    // can only be one adjacent facet at edge
  DLIList<CubitFacet*> pt_facets;
  point1->shared_facets(point2, pt_facets);
  if ( pt_facets.size() != 2 || !pt_facets.move_to(this) )
    return CUBIT_FAILURE;
  
    // get other facet
  CubitFacetData* other_facet = NULL;
  if( pt_facets.get() == this )
  {
    other_facet = dynamic_cast<CubitFacetData*>(pt_facets.next());
  }
  else
    if( pt_facets.next() == this )
    {
      other_facet = dynamic_cast<CubitFacetData*>( pt_facets.get() );
    }
    else
    {
      assert(0);
      return CUBIT_FAILURE;
    }

  assert(other_facet);
  
    // get indices on other facet
  int other_pt1_index = other_facet->point_index(point1);
  int other_pt2_index = (other_pt1_index+1)%3;
  int other_edge_index = (other_pt1_index+2)%3;
  if ( other_facet->point(other_pt2_index) != point2 ) {
    other_pt2_index = other_edge_index;
    other_edge_index = (other_pt1_index+1)%3;
  }
  assert( other_facet->point(other_pt2_index) == point2 );
  
    // check facet orientation
  int this_flip_use = this->edge_use(this_edge_index);
  int other_flip_use = other_facet->edge_use(other_edge_index);
  if ( this_flip_use == other_flip_use )
  {
    // Facets don't have consistant normals!!
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // get the opposite points on facets
  CubitPoint* this_other_pt = this->point(this_edge_index);
  CubitPoint* other_other_pt = other_facet->point(other_edge_index);
  if(this_other_pt == other_other_pt){
    PRINT_WARNING("Unable to perform flip.\n");
    return CUBIT_FAILURE;
  }
  
  
    // get the edge that is to be moved from this to the other facet
  CubitFacetEdge* this_trade_edge = this->edge(this_pt2_index);
    // get the edge thatis to be moved from the other facet to this
  CubitFacetEdge* other_trade_edge = other_facet->edge(other_pt1_index);
  if(this_trade_edge == other_trade_edge){
    PRINT_WARNING("Unable to perform flip (2).\n");
    return CUBIT_FAILURE;
  }
  int this_trade_use = this->edge_use(this_pt2_index);
  if ( this_trade_edge )
  {
    this_trade_edge->remove_facet(this);
    this_trade_edge->add_facet(other_facet);
  }
  
  int other_trade_use = other_facet->edge_use(other_pt1_index);
  if ( other_trade_edge )
  {
    other_trade_edge->remove_facet(other_facet);
    other_trade_edge->add_facet(this);
  }
  
    // get the edge to flip and change its points
  CubitFacetEdgeData* flip_edge 
    = dynamic_cast<CubitFacetEdgeData*>(edge(this_edge_index));
  if ( flip_edge )
  {
      // orient edge such that the edge uses stay the same
    int dir = (flip_edge->point(0) == point1);
    flip_edge->set_point( this_other_pt, 1-dir );
    flip_edge->set_point( other_other_pt, dir );
  }
  

    // change this facet
  this->edge( other_trade_edge, this_edge_index );
  this->edge_use( other_trade_use, this_edge_index );
  this->edge( flip_edge, this_pt2_index );
  this->edge_use( this_flip_use, this_pt2_index );
  point1->remove_facet(this);
  other_other_pt->add_facet(this);
  this->set_point( other_other_pt, this_pt1_index );
  
  
    // change the other facet
  other_facet->edge( this_trade_edge, other_edge_index );
  other_facet->edge_use( this_trade_use, other_edge_index );
  other_facet->edge( flip_edge, other_pt1_index );
  other_facet->edge_use( other_flip_use, other_pt1_index );
  point2->remove_facet(other_facet);
  this_other_pt->add_facet(other_facet);
  other_facet->set_point( this_other_pt, other_pt2_index );
  
  // make sure everything is correct
#ifndef NDEBUG
  for ( int i = 0; i < 3; i++ )
  {
    if ( this->edge(i) )
    {
      int start_index, end_index;
      if ( this->edge_use(i) == 1 )
      {
        start_index = (i+1)%3;
        end_index = (i+2)%3;
      }
      else
      {
        assert(this->edge_use(i) == -1);
        start_index = (i+2)%3;
        end_index = (i+1)%3;
      }
      assert ( this->edge(i)->point(0) == this->point(start_index) );
      assert ( this->edge(i)->point(1) == this->point(end_index) );
    }
    
    if ( other_facet->edge(i) )
    {
      int start_index, end_index;
      if ( other_facet->edge_use(i) == 1 )
      {
        start_index = (i+1)%3;
        end_index = (i+2)%3;
      }
      else
      {
        assert(other_facet->edge_use(i) == -1);
        start_index = (i+2)%3;
        end_index = (i+1)%3;
      }
      assert ( other_facet->edge(i)->point(0) == other_facet->point(start_index) );
      assert ( other_facet->edge(i)->point(1) == other_facet->point(end_index) );
    }
  }
#endif
  
  return CUBIT_SUCCESS;
}
