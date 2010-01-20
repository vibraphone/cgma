//-------------------------------------------------------------------------
// Filename      : Faceter.cpp
//
// Purpose       : Implementation of Faceter
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 03/01/02
//-------------------------------------------------------------------------

#include "Faceter.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "FaceterPointData.hpp"
#include "FaceterFacetData.hpp"
#include "CubitVector.hpp"

#include "DLIList.hpp"

#include "GeometryDefines.h"
#include "CubitDefines.h"
#include "GfxDebug.hpp"
#include "GeometryQueryEngine.hpp"
#include "GMem.hpp"

#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "CoEdge.hpp"
#include "RefVertex.hpp"
#include "CubitPlane.hpp"
#include "Curve.hpp"
#include "PointGridSearch.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 03/01/02
//-------------------------------------------------------------------------
Faceter::Faceter( RefFace *face ) 
    : thisRefFacePtr(face), gridCellScale(2.5)
{
  globalPointList = new DLIList <CubitPoint*>;
  gridSearchPtr = NULL;
  avoidedOverlap = CUBIT_FALSE;
}
//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 03/01/02
//-------------------------------------------------------------------------
Faceter::~Faceter()
{
  if (gridSearchPtr)
    delete gridSearchPtr;
  if ( globalPointList )
    delete globalPointList;
}
//-------------------------------------------------------------------------
// Purpose       : facet_surface: facets the surface.
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 03/01/02
//-------------------------------------------------------------------------
CubitStatus Faceter::facet_surface(DLIList <CubitFacet*> &results,
                                   DLIList <CubitPoint*> &point_list)
{
  if ( DEBUG_FLAG(129) )
  {
    GfxDebug::clear();
    GfxDebug::draw_ref_face_edges(thisRefFacePtr);
    GfxDebug::flush();
    int debug = 0;
    if ( debug )
    {
      GfxDebug::mouse_xforms();
      GfxDebug::flush();
    }
  }
  if ( thisRefFacePtr->number_of_Loops() > 1 )
    return CUBIT_FAILURE;
    //Get the ordered boundary loops.
  int ii, jj;
  DLIList <DLIList<CubitPoint*>*>  boundary_point_loops;
  DLIList <CubitPoint*> *tmp_list_ptr;
  CubitStatus stat = get_boundary_points( boundary_point_loops );
  if ( stat != CUBIT_SUCCESS )
  {
      //clean up the data...
    for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
    {
      tmp_list_ptr = boundary_point_loops.pop();
      for ( jj = tmp_list_ptr->size(); jj > 0; jj-- )
        delete tmp_list_ptr->pop();
      delete tmp_list_ptr;
    }
    return stat;
  }
    //Set up the gridsearch.
  double ratio = gridCellScale, cell_size = 0.0;
  max_min_edge_ratio(boundary_point_loops, ratio, cell_size);
  if (ratio <= gridCellScale) {
    ratio = gridCellScale;
  }
    //Get all of the points into a single list.
  for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
  {
    tmp_list_ptr = boundary_point_loops.get_and_step();
    for ( jj = tmp_list_ptr->size(); jj > 0; jj-- )
    {
      globalPointList->append(tmp_list_ptr->get_and_step());
    }
  }
  gridSearchPtr = new PointGridSearch(*globalPointList,
                                      cell_size,
                                      ratio);
      //fill in the grid...
  for ( ii = globalPointList->size(); ii > 0; ii-- )
    gridSearchPtr->add_point(globalPointList->get_and_step());

    //Now start faceting.
  stat = facet_loop( boundary_point_loops.get(), results );
    //clean up the data...
  for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
    delete boundary_point_loops.pop();
  if ( stat != CUBIT_SUCCESS )
  {
      //clean the data and return..
    for ( ii = results.size(); ii > 0; ii-- )
      delete results.pop();
    for ( ii = globalPointList->size(); ii > 0; ii-- )
      delete globalPointList->pop();
    return stat;
  }
    //Didn't add any points...
  point_list += *globalPointList;
  return CUBIT_SUCCESS;
}


CubitStatus Faceter::facet_loop(DLIList <CubitPoint*> *loop_ptr,
                                DLIList <CubitFacet*> &results)
{
    //Assume there is just this loop and facet it as if there was nothing
    //in the middle.  This is about the most stupid, simplistic method
    //of faceting I could think of.  I'm sure anyone reading this could
    //think of something more sophisticated and better.

    //First get the inteior angle for all of the points.
  DLIList <FaceterPointData*> ordered_point_list;
  int ii;
  double angle;
  CubitStatus stat;
  CubitPoint *curr_point;
  FaceterPointData *curr_faceter_data, *prev, *next;
  loop_ptr->reset();
  int debug = 0;
  for ( ii = 0; ii < loop_ptr->size(); ii++ )
  {
    stat = interior_angle(loop_ptr, angle);
    if ( stat == CUBIT_FAILURE )
      return CUBIT_FAILURE;
    curr_point = loop_ptr->get_and_step();
    if ( debug )
    {
      PRINT_INFO("%d\n",curr_point->id());
      GfxDebug::draw_point( curr_point->coordinates(), CUBIT_ORANGE);
    }
    curr_faceter_data = (FaceterPointData*) curr_point;
    prev = (FaceterPointData*) loop_ptr->prev(2);
    next = (FaceterPointData*) loop_ptr->get();

    if ( angle < CUBIT_RESABS || angle > DEGREES_TO_RADIANS( 359.95 ))
    {
      
      PRINT_ERROR("Zero degree angle in loop\n");
      RefEdge *owner_e = (RefEdge*) curr_faceter_data->owner();
      RefEdge *prev_e = (RefEdge*) prev->owner();
      RefEdge *next_e = (RefEdge*) next->owner();
      PRINT_ERROR("On edge %d, with next %d and prev %d\n", owner_e->id(),
                  next_e->id(), prev_e->id());
      return CUBIT_FAILURE;
    }
    curr_faceter_data->set_interior_angle(angle);
    curr_faceter_data->set_prev(prev);
    curr_faceter_data->set_next(next);
    ordered_point_list.append(curr_faceter_data);
  }
  ordered_point_list.sort(FaceterPointData::sort_by_angle);

    //Now go through and just cut off the smallest triangles...
    //The list should be ordered largest angle to smallest.
    //so go from the end to the front.
  CubitFacet *new_facet = NULL;
  int safe_size = 4* ordered_point_list.size();
  int safe_count = 0;
  int redo_count = ordered_point_list.size();
  while (ordered_point_list.size() >= 3 &&
         safe_count < safe_size )
  {
    curr_faceter_data = ordered_point_list.pop();
    new_facet = NULL;
    stat = facet_factory(curr_faceter_data, new_facet, ordered_point_list);
    if ( new_facet != NULL )
      results.append(new_facet);
    if ( stat != CUBIT_SUCCESS )
      return stat;
    if ( safe_count == redo_count )
    {
      if ( avoidedOverlap )
      {
        CubitBoolean reorder = CUBIT_FALSE;
        for ( ii = ordered_point_list.size(); ii > 0; ii-- )
        {
          curr_faceter_data = ordered_point_list.get_and_step();
          angle = curr_faceter_data->get_interior_angle();
          if ( angle > 2*CUBIT_PI )
          {
            angle = angle - 2*CUBIT_PI;
            curr_faceter_data->set_interior_angle(angle);
            reorder = CUBIT_TRUE;
          }
        }
        if ( reorder )
          ordered_point_list.sort(FaceterPointData::sort_by_angle);
          //set the next redo to be later.
        avoidedOverlap = CUBIT_FALSE;
      }
      redo_count += redo_count;
    }
    safe_count++;
  }
  if ( ordered_point_list.size() < 3 )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}
CubitStatus Faceter::facet_factory(FaceterPointData *curr_faceter_data,
                                   CubitFacet *&new_facet,
                                   DLIList <FaceterPointData*> &order_list)
{
  double angle;
  DLIList<CubitPoint*> loop_list;
  new_facet = NULL;
  if ( avoid_facet(curr_faceter_data, order_list) )
    return CUBIT_SUCCESS;
  FaceterPointData *prev = curr_faceter_data->get_prev();
  FaceterPointData *next = curr_faceter_data->get_next();
    //create the new facet.
    //create a triangle from these three points.
  new_facet = (CubitFacet*) new FaceterFacetData((CubitPoint*)prev,
                                                 (CubitPoint*)curr_faceter_data,
                                                 (CubitPoint*)next);
  CubitVector facet_normal = new_facet->normal();
  CubitVector surface_normal = thisRefFacePtr->normal_at(curr_faceter_data->coordinates());
  double dot = facet_normal%surface_normal;
  if ( dot < 0 )
  {
      //new need to try and do this in a different way...
    delete new_facet;
    new_facet = NULL;
    if ( order_list.size() < 3 )
    {
      PRINT_DEBUG_129("Last triangle is inverted. (surface %d)\n",
                      thisRefFacePtr->id());
      return CUBIT_FAILURE;
    }
      //Try and create facets on either side of this facet.
    return CUBIT_FAILURE;
  }
  
  int debug = 1;
  if ( debug )
    GfxDebug::draw_facet(new_facet, CUBIT_RED);
    //Now fix up the ordered_list.
    //First update the angles.
  prev->set_next(next);
  next->set_prev(prev);
  loop_list.clean_out();
  loop_list.append(prev->get_prev());
  loop_list.append(prev);
  loop_list.append(next);
  loop_list.step();
  interior_angle(&loop_list, angle);
  prev->set_interior_angle(angle);
  loop_list.clean_out();
  loop_list.append(prev);
  loop_list.append(next);
  loop_list.append(next->get_next());
  loop_list.step();
  interior_angle(&loop_list, angle);
  next->set_interior_angle(angle);
    //Now do basically a sort again.  There is probably
    //a faster way to do this but there are also much slower
    //ways in worst case too. At least this has guaranteed O(nlogn).
  order_list.sort(FaceterPointData::sort_by_angle);
  CubitPoint *tmp_point = (CubitPoint*) curr_faceter_data;
  gridSearchPtr->remove_point(tmp_point);
  return CUBIT_SUCCESS;
}
CubitBoolean Faceter::avoid_facet(FaceterPointData *curr_faceter_data,
                                  DLIList<FaceterPointData*> &order_list)
{
  
  int ii;
  CubitPoint *tmp_point;
  FaceterPointData *prev = curr_faceter_data->get_prev();
  FaceterPointData *next = curr_faceter_data->get_next();
    //get the closest edges to this point.
  CubitVector curr_v = curr_faceter_data->coordinates();
  CubitVector prev_v = prev->coordinates();
  CubitVector next_v = next->coordinates();
  double r1 = (curr_v - prev_v).length();
  double r2 = (curr_v - next_v).length();
  double radius;
  double angle;
  if ( r1 > r2 )
    radius = r1;
  else
    radius = r2;
  gridSearchPtr->set_neighborhood_bounds(curr_v, radius);
  DLIList<CubitPoint*> neighborhood_points;
  DLIList<CubitPoint*> close_points, loop_list;
  gridSearchPtr->get_neighborhood_points(neighborhood_points);
    //Loop through and see if there are nodes in here closer
    //than the r1 or r2.
  FaceterPointData *tmp_faceter;
  CubitVector tmp_v;
  radius *= radius;
  double angle1, angle2;
  for (ii = neighborhood_points.size(); ii > 0; ii-- )
  {
    tmp_point = neighborhood_points.get_and_step();
    tmp_faceter = (FaceterPointData*) tmp_point;
    if ( tmp_faceter == curr_faceter_data ||
         tmp_faceter == prev ||
         tmp_faceter == next ||
         tmp_faceter == next->get_next() ||
         tmp_faceter == prev->get_prev() )
      continue;
    tmp_v = tmp_faceter->coordinates();
    if ( (curr_v-tmp_v).length_squared() < radius )
      close_points.append(tmp_point);
  }
    //Now of these close points if any of them make
    //good trinagles for the other nodes then don't do
    //this triangle...
  CubitBoolean abort_facet = CUBIT_FALSE;
  for ( ii = close_points.size(); ii > 0; ii-- )
  {
    tmp_point = close_points.get_and_step();
	loop_list.clean_out();
    loop_list.append(tmp_point);
    loop_list.append((CubitPoint*)prev);
    loop_list.append((CubitPoint*)curr_faceter_data);
    loop_list.step();
    interior_angle(&loop_list, angle1);
    loop_list.clean_out();
    loop_list.append((CubitPoint*)curr_faceter_data);
    loop_list.append((CubitPoint*)next);
    loop_list.append(tmp_point);
    loop_list.step();
    interior_angle(&loop_list, angle2);
    if ( angle1 < DEGREES_TO_RADIANS(160) &&
         angle2 < DEGREES_TO_RADIANS(160))
    {
      int tmp_debug = 1;
      if ( tmp_debug )
      {
        GfxDebug::draw_point(prev->coordinates(), CUBIT_GREEN);
        GfxDebug::draw_point(curr_faceter_data->coordinates(), CUBIT_RED);
        GfxDebug::draw_point(tmp_point->coordinates(), CUBIT_BLUE);
        GfxDebug::draw_point(next->coordinates(), CUBIT_MAGENTA);
        loop_list.clean_out();
        loop_list.append(tmp_point);
        loop_list.append((CubitPoint*)prev);
        loop_list.append((CubitPoint*)curr_faceter_data);
        loop_list.step();
		interior_angle(&loop_list, angle1);
      }
      abort_facet = CUBIT_TRUE;
      break;
    }
  }
  if ( abort_facet )
  {
      //We need to add curr_facerter_data back to the list.
    order_list.insert(curr_faceter_data);
      //Set the angle above 360...
    angle = curr_faceter_data->get_interior_angle();
    curr_faceter_data->set_interior_angle(angle + 2*CUBIT_PI);
    avoidedOverlap = CUBIT_TRUE;
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitStatus Faceter::interior_angle(DLIList <CubitPoint*> *loop_ptr,
                                    double &my_angle)
{
  CubitPoint *current_node = loop_ptr->get();

  CubitVector point = loop_ptr->get()->coordinates();
  CubitVector to_prev = loop_ptr->prev()->coordinates();
  to_prev  -= point;
  CubitVector to_next = loop_ptr->next()->coordinates();
  to_next -= point;
  CubitVector surf_norm = thisRefFacePtr->normal_at ( point );
  if ( surf_norm.length_squared() < CUBIT_RESABS )
  {
      //Try getting it at one of the other nodes...
    surf_norm = thisRefFacePtr->normal_at(loop_ptr->next()->coordinates() );
    if (surf_norm.length_squared() < CUBIT_RESABS )
    {
        //Try getting it at one of the other nodes...
      surf_norm = thisRefFacePtr->normal_at(loop_ptr->prev()->coordinates() );
      if (surf_norm.length_squared() < CUBIT_RESABS )
      {
        PRINT_ERROR("Trying to get normal at point %d on surf %d.\n"
                    "       Normal length being returned equals zero.\n",
                    loop_ptr->get()->id(), thisRefFacePtr->id() );
        return CUBIT_FAILURE;
      }
    }
  }
  double angle = surf_norm.vector_angle ( to_next, to_prev );

    //Now if the angle is very small (zero) we don't know if this
    //should be 0 or 2PI, so do some tests.
  const double NEAR_ZERO = DEGREES_TO_RADIANS( 0.15 );
  const double NEAR_360  = DEGREES_TO_RADIANS( 359.95 );
  if( NEAR_ZERO < angle && angle < NEAR_360 )
  {
    my_angle = angle;
    return CUBIT_SUCCESS;
  }
    //If we are close to zero or 2PI then take some other loop nodes
    //into consideration to see if we can get a better approximation
    //as to we are zero or 2 PI.
  current_node = loop_ptr->get();
  CubitPoint *prev_node = loop_ptr->prev();
  CubitPoint *next_node = loop_ptr->next();
  int stop_it = loop_ptr->size();
  int count = 0;
  while (( angle < NEAR_ZERO || angle > NEAR_360 ) && ( count < stop_it ) )
  {
    prev_node = loop_ptr->prev( count + 2 );
    next_node = loop_ptr->next( count + 2 );
    to_prev = prev_node->coordinates() - point;
    to_next = next_node->coordinates() - point;
    angle = surf_norm.vector_angle ( to_next, to_prev );
    count++;
  }
  const double TWO_PI = 2.0 * CUBIT_PI;
  if ( angle < NEAR_ZERO )
    my_angle = TWO_PI;
  else if( angle < CUBIT_PI )
    my_angle = 0.0;
  else
    my_angle = TWO_PI;
  return CUBIT_SUCCESS;
}
void Faceter::max_min_edge_ratio( DLIList <DLIList <CubitPoint*>*> &boundary_point_loops,
                                  double &ratio,
                                  double &cell_size)
{
  double max_edge_length = CUBIT_DBL_MIN;
  double min_edge_length = CUBIT_DBL_MAX;
  double sum = 0.0;
  int total_count = 0;
  DLIList<CubitPoint*> *loopPtr;
  CubitPoint *node1, *node2;
  CubitVector edge;
  for (int i = boundary_point_loops.size(); i > 0; i--) {
    loopPtr = boundary_point_loops.get_and_step();
    node2 = loopPtr->prev();
    for (int j = loopPtr->size(); j > 0; j--) {
      node1 = node2;
      node2 = loopPtr->get_and_step();
      CubitVector p1 = node1->coordinates();
      CubitVector p2 = node2->coordinates();
      edge = p2-p1;
      double edge_length = edge.length_squared();
      if (edge_length > max_edge_length) max_edge_length = edge_length;
      if (edge_length < min_edge_length) min_edge_length = edge_length;
      total_count++;
      sum += edge_length;
    }
  }

  if (min_edge_length > CUBIT_RESABS) {
    ratio = sqrt(max_edge_length/min_edge_length);
  }
  else {
    ratio = sqrt(max_edge_length);
  }
  if ( total_count > 0 && sum > 0 )
    cell_size = sqrt(sum/total_count);
  else
    cell_size = 0.0;
}


CubitStatus Faceter::get_boundary_points( DLIList <DLIList<CubitPoint*>*> &boundary_point_loops ) const
{
  DLIList<DLIList<CoEdge*>*> co_edge_loops;
  thisRefFacePtr->co_edge_loops(co_edge_loops);
  int ii, jj;
  DLIList <CoEdge*> *co_edge_list_ptr;
  DLIList <CubitPoint*> *new_point_loop_ptr, tmp_point_list;
  RefEdge *ref_edge_ptr;
  CoEdge *co_edge_ptr;
  CubitStatus stat;
  CubitSense sense;
  
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list_ptr = co_edge_loops.get_and_step();
    new_point_loop_ptr = new DLIList <CubitPoint*>;
    for ( jj = co_edge_list_ptr->size(); jj > 0; jj-- )
    {
      co_edge_ptr = co_edge_list_ptr->get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      tmp_point_list.clean_out();
      stat = get_curve_facets( ref_edge_ptr, tmp_point_list );
      PRINT_DEBUG_129("curve %d has %d points\n",
                      ref_edge_ptr->id(),
                      tmp_point_list.size());
      if ( !stat )
        return CUBIT_FAILURE;
      tmp_point_list.reset();
        //the points are in order from start vertex to end vertex.
        //append them now according to the loop.
      sense = co_edge_ptr->get_sense();
      if ( CUBIT_FORWARD != sense )
        tmp_point_list.reverse();
        //Now take off the last point as it is a duplicate with the
        //other list...
      tmp_point_list.reset();
      delete tmp_point_list.pop();
      (*new_point_loop_ptr) += tmp_point_list;
    }
    CubitVector curr, prev;
    for ( jj = new_point_loop_ptr->size(); jj > 0; jj-- )
    {
      prev = new_point_loop_ptr->prev()->coordinates();
      curr = new_point_loop_ptr->get_and_step()->coordinates();
      if ( prev.about_equal(curr) )
      {
        PRINT_DEBUG_129("Points within tolerance in boundaryloop.\n");
        new_point_loop_ptr->back();
        delete new_point_loop_ptr->remove();
      }
    }
    boundary_point_loops.append(new_point_loop_ptr);
  }
  
    //clean up the list memory.
  for(ii = co_edge_loops.size(); ii>0; ii-- )
    delete co_edge_loops.pop();
  co_edge_loops.clean_out();
      
  return CUBIT_SUCCESS;
}

CubitStatus Faceter::get_curve_facets( RefEdge* curve, DLIList<CubitPoint*>& segments ) const
{
//const double COS_ANGLE_TOL =  0.965925826289068312213715; // cos(15)
  const double COS_ANGLE_TOL =  0.984807753012208020315654; // cos(10)
//const double COS_ANGLE_TOL =  0.996194698091745545198705; // cos(5)
  GMem curve_graphics;
  const double dist_tol = GEOMETRY_RESABS;
  const double dist_tol_sqr = dist_tol*dist_tol;
  Curve* curve_ptr = curve->get_curve_ptr();
  curve_ptr->get_geometry_query_engine()->get_graphics( curve_ptr, &curve_graphics );
  
  GPoint* gp = curve_graphics.point_list();
  CubitPoint* last = (CubitPoint*) new FaceterPointData( gp->x, gp->y, gp->z );
  ((FaceterPointData*)last)->owner(dynamic_cast<RefEntity*>(curve));
  CubitVector lastv = last->coordinates();
  segments.append( last );
  GPoint* end = gp + curve_graphics.pointListCount - 1;
  
  for( gp++; gp < end; gp++ )
  {
    CubitVector pos(  gp->x, gp->y, gp->z );
    CubitVector step1 = (pos - lastv);
    double len1 = step1.length();
    if( len1 < dist_tol ) continue;
    
    GPoint* np = gp + 1;
    CubitVector next( np->x, np->y, np->z );
    CubitVector step2 = next - pos;
    double len2 = step2.length();
    if( len2 < dist_tol ) continue;
    
    double cosine = (step1 % step2) / (len1 * len2);
    if( cosine > COS_ANGLE_TOL ) continue;
    
    last = new FaceterPointData( pos );
    ((FaceterPointData*)last)->owner(dynamic_cast<RefEntity*>(curve));
    segments.append( last );
    lastv = last->coordinates();
  }
  
  CubitVector last_pos( gp->x, gp->y, gp->z );
  segments.last();
  while( (last_pos - (segments.get()->coordinates())).length_squared() < dist_tol_sqr )
  {
    delete segments.pop();
    segments.last();
  }
  CubitPoint *tmp_point = (CubitPoint*) new FaceterPointData( last_pos );
  segments.append( tmp_point );
  ((FaceterPointData*)tmp_point)->owner( dynamic_cast<RefEntity*>(curve) );
    
  // Now check if the segment list is reversed wrt the curve direction.
  segments.reset();
  double u1, u2;
  if( segments.size() > 2 )
  {
    u1 = curve->u_from_position( (segments.next(1)->coordinates()) );
    u2 = curve->u_from_position( (segments.next(2)->coordinates()) );    
  }
  else
  {
    u1 = curve->u_from_position( (segments.get()->coordinates() ) );
    u2 = curve->u_from_position( (segments.next()->coordinates()) );
  }
  if( (u2 < u1) && (curve->start_param() <= curve->end_param()) )
    segments.reverse();

    //Make sure we don't have duplicate points.
  int jj;
  CubitVector curr, prev;
  for ( jj = segments.size(); jj > 0; jj-- )
  {
    prev = segments.prev()->coordinates();
    curr = segments.get_and_step()->coordinates();
    if ( prev.about_equal(curr) )
    {
      PRINT_DEBUG_129("Points on curve %d within tolerance...\n", curve->id());
      segments.back();
      delete segments.remove();
    }
  }
  return CUBIT_SUCCESS;
}

