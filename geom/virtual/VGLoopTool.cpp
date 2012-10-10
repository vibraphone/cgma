
#include "VGLoopTool.hpp"
#include "GMem.hpp"
#include "GeometryQueryEngine.hpp"
#include "CubitVector.hpp"
/*
//-------------------------------------------------------------------------
// Purpose       : Remove a curve from a loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
remove_curve( COEDGE* coedge1, COEDGE* coedge2,
              DLIList<COEDGE*>& new_loop_coedges )
{
  new_loop_coedges.clean_out();
  
  if( !coedge1->get_loop() || !coedge2->get_loop() )
    return CUBIT_FAILURE;
  
    // remove sipe or split loop?
  if( coedge1->get_loop() == coedge2->get_loop() )
  {
      // split loop?
    if( coedge1->next() != coedge2 && coedge1->prev() != coedge2 )
    {
      COEDGE* coedge = coedge1->next();
      while( coedge != coedge2 )
      {
        COEDGE* next = coedge;
        coedge->get_loop()->remove( coedge );
        new_loop_coedges.append( coedge );
        coedge = next;
      }
    }
    
    coedge1->get_loop()->remove( coedge1 );
    coedge2->get_loop()->remove( coedge2 );
  }
  
    // stitch/combine loops
  else
  {
    COEDGE* coedge;
    
      // insert coedges
    while( coedge2->next() != coedge2 )  // loop has more than 1 coedge
    {
      coedge = coedge2->next();
      coedge2->get_loop()->remove( coedge );
      coedge1->insert_before( coedge, coedge1 );
    }
    coedge1->get_loop()->remove( coedge1 );
    coedge2->get_loop()->remove( coedge2 );
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Insert a curve in a loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
insert_curve( COEDGE* coedge1, COEDGE* coedge2,
              SURFACE* surface, LOOP*& new_loop1, LOOP*& new_loop2 )
{
  assert( coedge1 && coedge2 && 
          coedge1->get_curve() == coedge2->get_curve() &&
          coedge1->sense() != coedge2->sense() );
  new_loop1 = new_loop2 = 0;
  
  if( coedge1->sense() == CUBIT_REVERSED )
  {
    COEDGE* tmp = coedge1;
    coedge2 = coedge2;
    coedge2 = tmp;
  }
  
    // Find loop insertion locations.
  COEDGE *start_prev = 0, *end_prev = 0;
  if( ! previous_coedge( coedge2, surface, start_prev ) ||
      ! previous_coedge( coedge1, surface, end_prev ) )
    return CUBIT_FAILURE;
  
    // If no loop was passed, just pass back the coedges
    // so that the caller can create new loop(s)
  if( !(start_prev || end_prev) )
  {
    new_loop1 = new LOOP;
    new_loop1->insert_after( coedge1, 0 );
    if( coedge1->start_point() == coedge1->end_point() )
    {
      new_loop2 = new LOOP;
      new_loop2->insert_after( coedge2, 0 );
    }
    else
    {
      new_loop1->insert_after( coedge2, coedge1 );
    }
    
    return CUBIT_SUCCESS;
  }
  
  
    // sipe?
  if( !loop1_prev || !loop2_prev )
  {
    COEDGE *prev = loop1_prev ? loop1_prev : loop2_prev;
    if( coedge1->start_point() == prev->end_point() )
    {
      prev->get_loop()->insert_after( coedge1, prev );
      prev->get_loop()->insert_after( coedge2, coedge1 );
    }
    else if( coedge1->end_point() == prev->end_point() )
    {
      prev->get_loop()->insert_after( coedge2, prev );
      prev->get_loop()->insert_after( coedge1, coedge2 );
    }
    else
    {
      assert(0);
      return CUBIT_FAILURE;
    }
  }
  
    // combine loops?
  else if( loop1_prev->get_loop() != loop2_prev->get_loop() )
  {
    COEDGE* prev = 0;

      // Which of the two coedges for the curve are
      // inserting do we want to begin with?
    if( coedge1->start_point() == loop1_prev->end_point() )
    {
      loop1_prev->get_loop()->insert_after( coedge1, loop1_prev );
      loop1_prev->get_loop()->insert_after( coedge2, coedge1 );
      prev = coedge1;
    }
    else if( coedge2->start_point() == loop1_prev->end_point() )
    {
      loop1_prev->get_loop()->insert_after( coedge2, loop1_prev );
      loop1_prev->get_loop()->insert_after( coedge1, coedge2 );
      prev = coedge2;
    }
    else
    {
      assert(0);
      return CUBIT_FAILURE;
    }
    
    COEDGE* coedge = loop2_prev->next();
    COEDGE* next = 0;
    while( loop2_prev->get_loop() != loop1_prev->get_loop() )
    {
      next = coedge->next();
      coedge->get_loop()->remove( coedge );
      loop1_prev->get_loop()->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
  }
  
    // split the loop
  else
  {
    assert( loop1_prev->get_loop() == loop2_prev->get_loop() );
    
      // check for and handle a hole insersecting the original
      // loop at just one vertex
    if( coedge1->start_point() == coedge1->end_point() )
    {
      assert( loop1_prev == loop2_prev );
      
      CubitVector prev_tan, coe1_tan, coe2_tan, norm, junk;
      CubitVector point = coedge1->start_point()->coordinates();
      loop1_prev->get_curve()->closest_point( point, junk, &prev_tan );
      if( loop1_prev->sense() == CUBIT_FORWARD ) // yes, I mean forward here!!
        prev_tan *= -1.0;
      coedge1->get_curve()->closest_point( point, junk, &coe1_tan );
      coe2_tan = coe1_tan;
      if( coedge1->sense() == CUBIT_REVERSED )
        coe1_tan *= -1.0;
      if( coedge2->sense() == CUBIT_REVERSED )
        coe2_tan *= -1.0;
      loop1_prev->get_loop()->get_surface()->closest_point( point, 0, &norm );
      
      double angle1 = norm.vector_angle( prev_tan, coe1_tan );
      double angle2 = norm.vector_angle( prev_tan, coe2_tan );
      
      if( angle2 < angle1 )
      {
        COEDGE* temp = coedge2;
        coedge2 = coedge1;
        coedge1 = temp;
      }
    }
    
    else if( coedge1->end_point() == loop2_prev->end_point() )
    {
      COEDGE* temp = coedge2;
      coedge2 = coedge1;
      coedge1 = temp;
    }
  
    
    new_loop1 = new LOOP;
    loop1_prev->insert_after( coedge1, loop2_prev );
    new_loop1->insert_after( coedge2, 0 );
    COEDGE* prev = coedge2;
    
    COEDGE* coedge = loop2_prev->next();
    COEDGE* stop = loop1_prev->next();
    while( coedge != stop )
    {
      COEDGE* next = coedge->next();
      coedge->get_loop()->remove( coedge );
      new_loop1->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
  }
  
  return CUBIT_SUCCESS;
}
*/
/*
template <class SURFACE, class LOOP, class COEDGE, class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
VGLoopTool::curves_in_loop( LOOP* loop, POINT* point, DLIList<CURVE*>& results )
{
  CURVE* curve = 0;
  while ( (curve = point->next_curve(curve)) )
  {
    COEDGE* coedge = 0;
    while ( (coedge = curve->next_coedge(coedge)) )
    {
      if (coedge->get_loop() == loop)
      {
        results.append(curve);
        break;
      }
    }
  }
  return CUBIT_SUCCESS;
}

template <class SURFACE, class LOOP, class COEDGE, class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
VGLoopTool::traverse_coedges( LOOP* loop, COEDGE* coedge, DLIList<COEDGE*>& results )
{  
    // Traverse coedges, appending them to the list until
    // the end point is reached.
  POINT* loop_point = coedge->start_point();
  results.append(coedge);
  DLIList<CURVE*> loop_curves;
  while ( coedge->end_point() != loop_point )
  {
    POINT* point = coedge->end_point();
    loop_curves.clean_out();
    curves_in_loop( loop, point, loop_curves );
    
    COEDGE *boundary_coedge = 0, *split_coedge = 0;
    while( loop_curves.size() )
    {
      CURVE* curve = loop_curves.pop();
      if ( curve == coedge->get_curve() )
        continue;
      
      COEDGE* curve_coedge = 0;
      while ( (curve_coedge = curve->next_coedge(curve_coedge)) )
      {
        if ((curve_coedge->get_loop() != loop_to_split()) ||
            (curve_coedge->start_point() != point))
          continue;
        
        if (curve_coedge->get_curve()->mark == 2)
        {
          if (split_coedge) 
            { assert(0); return CUBIT_FAILURE; }
          split_coedge = curve_coedge;
        }
        else if(curve_coedge->get_curve()->mark == 0)
        {
          if (boundary_coedge) 
            { assert(0); return CUBIT_FAILURE; }
          boundary_coedge = curve_coedge;
        }
      }
    }
    
    if (split_coedge)
      coedge = split_coedge;
    else if(boundary_coedge)
      coedge = boundary_coedge;
    else
      { assert(0); return CUBIT_FAILURE; }

    results.append(coedge);
  }
  
  return CUBIT_SUCCESS;
}


template <class SURFACE, class LOOP, class COEDGE, class CURVE, class POINT>
LOOP* VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
VGLoopTool::split_loop( LOOP* loop_to_split )
{
    // Make sure all coedge curve marks are cleared,
    // and count coedges (use the count to make sure
    // things are going okay later)
  int coedge_count = 0;
  COEDGE* coedge = loop_to_split->first_coedge();
  do {
    coedge->mark = 0;
    coedge->get_curve()->mark = 0;
    coedge_count++;
    coedge = loop_to_split->next_coedge(coedge);
  } while ( coedge != loop_to_split->first_coedge() );
  
    // Identify non-manifold curves, marking them either
    // with a 2 if they can be used to split the surface
    // (if they are part of a connected patch for which the
    // bounadary of that patch intersects the surface boundary
    // at all points) or a 3 if they are other non-manifold
    // curves.  This will get a bit tricky if there are
    // non-manifold curves hanging off of the patch of
    // split curves.
  
    // First for all non-manifold curves, if the curve has
    // a point that is not shared with any other curve, mark
    // it with a 3.  Otherwise mark it with a 2.
  DLIList<CURVE*> curve_stack, loop_curves, nonman_curves;
  coedge = loop_to_split->first_coedge();
  do {
    CURVE* curve = coedge->get_curve();
      // If we haven't done this curve yet and it is non-manifold
    if ( !curve->mark && curve->find_next(coedge) )
    {
      curve->mark = 2;
 
      for ( int i = 2; i--; )
      {
        POINT* point = i ? curve->start_point() : curve->end_point();
        loop_curves.clean_out();
        curves_in_loop( loop_to_split, curve->start_point(), loop_curves );
        if ( loop_curves.size() == 1 )
        {
          curve->mark = 3;
          curve_stack.append(curve);
          nonman_curves.append(curve);
          break;
        }
      }
    }
    
    coedge = loop_to_split->next_coedge(coedge);
  } while( coedge != loop_to_split->first_coedge() );
  
    // Now for each curve we marked with a three, traverse
    // and mark adjacent curve until we come to a point
    // connected to more that two curves.
  while( curve_stack.size() ) 
  {
    CURVE* curve = curve_stack.pop();
    
    for ( int i = 2; i--; )
    {
      POINT* point = i ? curve->start_point() : curve->end_point();
      loop_curves.clean_out();
      curves_in_loop( loop_to_split, curve->start_point(), loop_curves );
      if ( loop_curves.size() != 2 )
        continue;
      
      loop_curves.move_to(curve);
      CURVE* other = loop_curves.next();
      assert(other && other != curve );
      
      if ( other->mark != 3 )
      {
        other->mark = 3;
        curve_stack.append(other);
        nonman_curves.append(other);
      }
    }
  }
  
    // Find a split curve to begin at.
  coedge = loop_to_split->first_coedge();
  do
  {
    if ( coedge->get_curve()->mark == 2 )
      break;
    coedge = loop_to_split->next_coedge(coedge);
  } while( coedge != loop_to_split->first_coedge() );
  
  if ( coedge->get_curve()->mark != 2 )
    { assert(0); return 0; }
    
    // Get two coedges for split curve
  COEDGE *coedge1 = coedge, *coedge2 = 0;
  coedge = 0;
  while ( (coedge = coedge1->get_curve()->next_coedge(coedge)) )
    if (coedge != coedge1)
      { coedge2 = coedge; break; }
  assert( coedge2 && coedge2->sense() != coedge1->sense() );
  
    // Get coedges for each loop by traversing from 
    // split curve coedges
  DLIList<COEDGE*> old_loop_coedges, new_loop_coedges;
  if ( !traverse_loop(loop_to_split, coedge1, old_loop_coedges) )
    return 0;
  if ( !traverse_loop(loop_to_split, coedge2, new_loop_coedges) )
    return 0;
  
    // Make sure everything adds up
  if ( old_loop_coedges.size() + new_loop_coedges.size() 
       + 2*nonman_curves.size() != coedge_count)
    { assert(0); return 0; }
  
    // Leave the larger list on the existing loop, as this
    // is more likely to result in the expected behavior
  if (new_loop_coedges.size() > old_loop_coedges.size())
  {
    DLIList<COEDGE*> tmp_list(old_loop_coedges);
    old_loop_coedges = new_loop_coedges;
    new_loop_coedges = tmp_list;
  }
  
    // Remove all coedges from old loop (faster this way)
  while( loop_to_split->first_coedge() )
    loop_to_split->remove( loop_to_split->first_coedge() );
  
    // Put old_loop_coedges in loop_to_split
  COEDGE* next = 0;
  while ( old_loop_coedges.size() )
  {
    coedge = old_loop_coedges.pop();
    coedge->get_curve()->mark = 0;
    loop_to_split->insert_before( coedge, next );
    next = coedge;
  }
  
    // Put new_loop_coedges in new_loop
  LOOP* new_loop = new LOOP;
  next = 0;
  while ( new_loop_coedges.size() )
  {
    coedge = new_loop_coedges.pop();
    coedge->get_curve()->mark = 0;
    new_loop->insert_before( coedge, next );
    next = coedge;
  }
  
    // Now sort out other non-manifold curves
  
    // Clear marks
  for ( int j = nonman_curves.size(); j--; )
    nonman_curves->step_and_get()->mark = 0;
  
  if( nonman_curves.size() )
    insert_nonmanifold_curves( nonman_curves, loop_to_split, new_loop );
  return new_shell;
}
  
template <class SURFACE, class LOOP, class COEDGE, class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
VGLoopTool::insert_nonmanifold_curves( DLIList<CURVE*>& curve_stack,
                                       LOOP* loop1, LOOP* loop2 )
{
  return CUBIT_FAILURE;
}

*/

/*  
//-------------------------------------------------------------------------
// Purpose       : Check if a curve occurs an odd number of times in a loop
//
// Special Notes : Input is the coedge connecting the relevent curve & loop
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
bool VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
odd_coedge_count( COEDGE* coedge )
{
  COEDGE* curve_coedge = 0;
  int count = 0;
  while( curve_coedge = coedge->get_curve()->next_coedge( curve_coedge ) )
    if( curve_coedge->get_loop() == coedge->get_loop() )
      count++;
    
  return (bool)(count % 2);
}

//-------------------------------------------------------------------------
// Purpose       : Find the previous coedge in loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/28/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
previous_coedge( COEDGE* coedge, SURFACE* surface, COEDGE*& result )
{
  result = 0;
  
  POINT* point = coedge->end_point();
  CubitVector normal, tanget, closest;
  surface->closest_point( point->coordinates(), &closest, &normal );
  
  coedge->get_curve()->closest_point( point->coordinates(), closest, &tangent );
  if( coedge->sense() == CUBIT_REVERSED )
    tangent *= -1.0;
  
  CURVE* curve = 0;
  while( curve = point->next_curve( curve ) )
  {
    COEDGE* coedge = 0;
    while( coedge = curve->next_coedge( coedge ) )
    {
      if( coedge->end_point() != point )
        continue;
      
      if( !coedge->get_loop() || coedge->get_loop()->get_surface() != surface )
        continue;
      
      if( ! result )
      {
        result = coedge;
      }
      else if( result->get_loop() != coedge->get_loop() )
      {
        result = 0;
        return CUBIT_FAILURE;
      }
      else
      {
        CubitVector result_tan, coedge_tan, closest;
        result->get_curve()->closest_point( point->coordinates(), closest, &result_tan );
        coedge->get_curve()->closest_point( point->coordinates(), closest, &coedge_tan );
        if( result->sense() == CUBIT_FORWARD )
          result_tan *= -1.0;
        if( coedge->sense() == CUBIT_FORWARD )
          result_tan *= -1.0;
      
        double result_angle = normal.vector_angle( result_tan, tangent );
        double coedge_angle = normal.vector_angle( coedge_tan, tangent );
        if( result_angle < 0.0 ) result_angle += 2.0*CUBIT_PI;
        if( coedge_angle < 0.0 ) coedge_angle += 2.0*CUBIT_PI;
        if( coedge_angle < result_angle )
          result = coedge;
      }
    }
  }
  
  return CUBIT_SUCCESS;
}

*/


//-------------------------------------------------------------------------
// Purpose       : Get a polygon representation of a loop as a point list
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
get_loop_polygon( COEDGE* first_coedge, DLIList<CubitVector*>& polygon )
{
  if( !first_coedge )
    return CUBIT_FAILURE;
  
  DLIList<CubitVector*> interior;
  CubitSense sense;
  polygon.clean_out();

  COEDGE* coedge = first_coedge;
  do
  {
      polygon.append( new CubitVector( coedge->start_point()->coordinates() ) );

      interior.clean_out();
      if( !coedge->get_curve()->get_interior_extrema( interior, sense ) )
      {
        while( polygon.size() ) delete polygon.pop();
        return CUBIT_FAILURE;
      }

      if( sense != coedge->sense() )
        interior.reverse();
      polygon += interior;
       
    coedge = coedge->next();
  } while( coedge != first_coedge );
  
  return CUBIT_SUCCESS;
}
  
  
//-------------------------------------------------------------------------
// Purpose       : Calculate loop angle metric
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
double VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
loop_angle_metric( COEDGE* first_coedge )
{
  if( !first_coedge )
    return 0.0;
    
  assert( first_coedge 
       && first_coedge->get_loop()
       && first_coedge->get_loop()->get_surface() );
  
  DLIList<CubitVector*> polygon;
  get_loop_polygon( first_coedge, polygon );
  
  double result = 0.0;
  CubitVector *point[3], segment[2], normal;
  point[0] = polygon.step_and_get();
  point[1] = polygon.step_and_get();
  segment[0] = *point[1] - *point[0];
  for( int i = polygon.size(); i--; )
  {
    point[2] = polygon.step_and_get();
    segment[1] = *point[1] - *point[2];
    first_coedge->get_loop()->get_surface()->closest_point( *point[1], 0, &normal );
    result += normal.vector_angle( segment[1], segment[0] );
    
    point[1] = point[2];
    segment[0] = -segment[1];
  }
  
  result /= CUBIT_PI;
  result -= polygon.size();
  
  while( polygon.size() )
    delete polygon.pop();
  
  return result;
}

template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
void VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
get_loop_polyline( COEDGE* first_coedge, std::vector<CubitVector>& list )
{
  COEDGE* coedge = first_coedge;
  GMem gmem;
  list.resize(0);
  
  do {
    coedge->get_curve()->get_geometry_query_engine()->
      get_graphics( coedge->get_curve(), &gmem );
    int num_pts = gmem.pointListCount;
    if ( num_pts > 1 ) {
    
      int start = list.size();
      list.resize( start + num_pts - 1 );
      std::vector<CubitVector>::iterator itor = list.begin() + start;
      if ( coedge->sense() == CUBIT_FORWARD ) {
        GPoint* pitor = gmem.point_list();
        GPoint* pend = pitor + num_pts - 1;
        for ( ; pitor != pend; ++pitor ) {
          itor->set( pitor->x, pitor->y, pitor->z );
          ++itor;
        }
      } else {
        GPoint* pend = gmem.point_list();
        GPoint* pitor = pend + num_pts - 1;
        for ( ; pitor != pend; --pitor ) {
          itor->set( pitor->x, pitor->y, pitor->z );
          ++itor;
        }
      }
   }
    
    coedge = coedge->next();
  } while ( coedge != first_coedge );
}
      

/*
//-------------------------------------------------------------------------
// Purpose       : Check if a position is to the "left" of a coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
bool VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
inside_of_curve( const CubitVector& tangent,
                 const CubitVector& curve_pos,
                 const CubitVector& surf_pos,
                 const CubitVector& normal )
{
  CubitVector cross = tangent * ( surf_pos - curve_pos );
  return cross % normal > -CUBIT_RESABS;
}

//-------------------------------------------------------------------------
// Purpose       : Find the coedge in a loop that is closest to the passed 
//                 point
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
COEDGE* VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
closest_loop_coedge( const CubitVector& from_pt,
                     COEDGE* first_coedge,
                     COEDGE*& other_coedge,
                     CubitVector* closest )
{
  CubitVector closest_pt, current_pt;
  double current_dst, closest_dst;
  COEDGE* closest_coedge = 0;
  other_coedge = 0;

    // Skip coedges for curves that occur twice in the loop (sipes, etc)
  COEDGE* coedge = first_coedge;
  while( !odd_coedge_count(coedge) )
  {
    coedge = coedge->next();
    if( coedge == first_coedge )
      return false;
  }

    // Search for a) a curve who's bounding box contains the point
    // or b) the curve who's box is closest to the point
  closest_coedge = coedge;
  closest_dst = coedge->get_curve()->bounding_box().distance_squared(from_pt);
  coedge = coedge->next();
  while( (coedge != first_coedge) && (closest_dst > CUBIT_RESABS) )
  {
    if( odd_coedge_count(coedge) )
    {
      current_dst = coedge->get_curve()->bounding_box().distance_squared(from_pt);
      if( current_dst < closest_dst )
      {
        closest_coedge = coedge;
        closest_dst = current_dst;
      }
    }
    coedge = coedge->next();
  }
  
    // Now find actual distance from the curve 
  closest_coedge->get_curve()->closest_point_trimmed( from_pt, closest_pt );
  closest_dst = (from_pt - closest_pt).length_squared();
  
    // Look for any curves that are closer
  first_coedge = closest_coedge;
  coedge = first_coedge.next();
  while( coedge != first_coedge )
  {
    if( odd_coedge_count(coedge) )
    {
      if( coedge->get_curve()->bounding_box().distance_squared(from_pt) < closest_dst )
      {
        coedge->get_curve()->closest_point_trimmed( from_pt, current_pt );
        current_dst = (from_pt - current_pt).length_squared();

        if( (closest_pt - current_pt).length_squared() < RESSQR )
        {
          if( current_dst < closest_dst )
          {
            closest_dst = current_dst;
            closest_pt = current_pt;
          }
          other_coedge = coedge;
        }
        else if( current_dst < closest_dst )
        {
          closest_coedge = coedge;
          closest_dst = current_dst;
          closest_pt = current_pt;
          other_coedge = 0;
        }
      }
    }
    coedge = coedge->next();
  }

  //make sure we have things in the correct order
  if( other_coedge )
  {
    if( closest_coedge->start_point() == other_coedge->end_point() )
    {
      COEDGE* tmp = closest_coedge;
      closest_coedge = other_coedge;
      other_coedge = tmp;
    }
    else if( closest_coedge->end_point() != other_coedge->start_point() )
    {
      other_coedge = 0;
    }
  }    

  if( closest ) *closest = closest_pt;
  return closest_coedge;
}
  

//-------------------------------------------------------------------------
// Purpose       : Test if a position is on the surface side of a loop.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
bool VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
is_position_within_loop( const CubitVector& position,
                         COEDGE* first_coedge,
                         CubitVector* closest_on_loop,
                         CubitVector* passed_normal )
{
  CubitVector surf_point, normal, pt_on_curve;
  if( passed_normal )
  {
    surf_point = position;
    normal = *pased_normal;
  }
  else
  {
    first_coedge->get_loop()->get_surface()
      ->closest_point( position, &surf_point, &normal );
  }
  
  COEDGE* other_coedge = 0;
  COEDGE* closest_coedge = 
    closest_loop_coedge( first_coedge, position, other_coedge, &pt_on_curve );
  if( !closest_coedge ) 
    return true;
  
  if( closest_on_loop )
    *closest_on_loop = pt_on_curve;
  
  CubitVector coe_normal, coe_cross, tangent1, tangent2, junk;
  
  if( ! other_coedge )
  {
    double u = closest_coedge->get_curve()->u_from_position( pt_on_curve );
    
    // Special case: closest point on loop at G1 discontinuity 
    if( closest_coedge->get_curve()->G1_discontinuous( u, &tangent1, &tangent2 ) )
    {
      if( closest_coedge->sense() == CUBIT_REVERSED )
      {
        tangent1 *= -1.0;
        tangent2 *= -1.0;
      }
      
      closest_coedge->get_loop()->get_surface()->
        closest_point( pt_on_curve, 0, &coe_normal );
      coe_cross = tangent1 * tangent2;
      bool inside1 = inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
      bool inside2 = inside_of_curve( tangent2, pt_on_curve, surf_point, normal );
      if( cross % coe_normal > -CUBIT_RESABS )
        return inside1 && inside2;
      else
        return inside1 || inside2;
    }
    
    else // **** This is the normal, non-special case part ****
    { 
      closest_coedge->closest_point( pt_on_curve, junk, &tangent1 );
      if( closest_coedge->sense() == CUBIT_REVERSED )
        tangent1 *= -1.0;
      
      return inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
    }
  }
    // Special case: closest point on loop at vertex 
  else
  {
    closest_coedge->get_curve()->closest_point( pt_on_curve, junk, &tangent1 );
      other_coedge->get_curve()->closest_point( pt_on_curve, junk, &tangent2 );
    if( closest_coedge->sense() == CUBIT_REVERSED )
      tangent1 *= -1.0;
    if(   other_coedge->sense() == CUBIT_REVERSED )
      tangent2 *= -1.0;
    
    closest_coedge->get_loop()->get_surface()->
      closest_point( pt_on_curve, 0, &coe_normal );
    
    coe_cross = tangent1 * tangent2;
    bool inside1 = inside_of_curve( tangent1, pt_on_curve, surf_point, normal );
    bool inside2 = inside_of_curve( tangent2, pt_on_curve, surf_point, normal );
    
    if( coe_cross % coe_normal > -CUBIT_RESABS )
      return inside1 && inside2;
    else
      return inside1 || inside2;
  }
}


//-------------------------------------------------------------------------
// Purpose       : Very rough approximation of the area bounded by a loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/14/02
//-------------------------------------------------------------------------
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
double VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
loop_area( COEDGE* first_coedge )
{
  int i;
  
  DLIList<CubitVector*> polygon;
  get_loop_polygon( first_coedge );

    // Get surface normal at each polygon point
  DLIList<CubitVector*> normals;
  polygon.reset();
  for( i = polygon.size(); i--; )
  {
    CubitVector norm, *pt = polygon.get_and_step();
    first_coedge->get_loop()->get_surface()->closest_point(*pt,0,&norm);
    normals.append( new CubitVector(norm) );
  }
  
    // Average surface normals
  normal_list.reset();
  CubitVector mean_normal(0.,0.,0.);
  for( i = normals.size(); i--; )
    mean_normal += *normals.get_and_step();
  mean_normal /= normals.size();
  
    // Choose base point at which the surface normal
    // is closest to the mean normal from above
  point_list.reset();
  normal_list.reset();
  CubitVector* base_point = polygon.get_and_step();
  CubitVector* base_normal = normals.get_and_step();
  double largest_dot = mean_normal % *base_normal;
  
  for( i = polygon.size(); i > 1; i-- )
  {
    CubitVector* normal = normals.get_and_step();
    double dot = mean_normal % *normal;
    if( dot > largest_dot )
    {
      base_point = polygon.get();
      base_normal = normal;
    }
    polygon.step();
  }
  
  polygon.move_to( base_point );
  normals.move_to( base_normal );
  polygon.step();
  normals.step();
  
    // Do simple 2-D area calculation and hope it isn't too
    // inaccurate for the 3-D refface.
  double double_area = 0.0;
  CubitVector *prev_point = polygon.get_and_step(),
              *curr_point = polygon.get_and_step();
  CubitVector *prev_normal = normals.get_and_step(),
              *curr_normal = normals.get_and_step();
    
  while( curr_point != base_point )
  {
    CubitVector cross = (*prev_point - *base_point)*(*curr_point  - *base_point);
    CubitVector normal = *base_normal + *prev_normal + *curr_normal;
    double dot_product = cross % normal;
    double_area += dot_product >= 0.0 ? cross.length() : -(cross.length());
    
    prev_point = curr_point;
    prev_normal = curr_normal;
    curr_point = polygon.get_and_step();
    curr_normal = normals.get_and_step();
  }

  while( polygon.size() > 0 ) delete polygon.pop();
  while( normals.size() > 0 ) delete normals.pop();
  return double_area / 2.0;
}
  
template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
closest_point_trimmed( DLIList<COEDGE*>& surface_coedges, const CubitVector& from_pt,
                       CubitVector& result_pt, CubitPointContainment* containment )
{
  CubitVector pt_on_surf, normal;
  if( ! surface->closest_point( from_pt, &pt_on_surf, &normal ) )
    return CUBIT_FAILURE;
  

  CubitVector pt_on_curve;  
  COEDGE *closest_coedge = 0, *other_coedge = 0;
  if( ! closest_surface_coedge( surface, from_pt, closest_coedge,
            other_coedge, &pt_on_curve ) )
    return CUBIT_FAILURE;
  
  CubitVector coe_normal, coe_cross, tangent1, tangent2, junk;
  bool inside;
  if( ! other_coedge )
  {
    CURVE* curve_ptr = closest_coedge->get_curve();
    double u = curve_ptr->u_from_position( pt_on_curve );
    
    // **** Special case: closest point on loop at G1 discontinuity ****
    if( curve_ptr->G1_discontinuous( u, &tangent1, &tangent2 ) )
    {
      if( closest_coedge->sense() == CUBIT_REVERSED )
      {
        tangent1 *= -1.0;
        tangent2 *= -1.0;
      }
      surface->closest_point( pt_on_curve, NULL, &coe_normal );
      coe_cross = tangent1 * tangent2;
      double sum  = (coe_cross + coe_normal).length_squared();
      double diff = (coe_cross - coe_normal).length_squared();
    
      CubitBoolean inside1 = 
        inside_of_curve( tangent1, pt_on_curve, pt_on_surf, normal );
      CubitBoolean inside2 = 
        inside_of_curve( tangent2, pt_on_curve, pt_on_surf, normal );
        
      if( (sum > diff) || (diff < CUBIT_DBL_MIN) ) 
        //discontinuity is at a convexity
      {
        //the point must be inside of both sub-edges
        inside = inside1 && inside2;
      }
      else //discontinuity is at a concavity
      {
        //the point must be inside of one of the sub-edges
        inside = inside1 || inside2;
      }
    }

    else // **** This is the normal, non-special case part ****
    { 
      curve_ptr->closest_point( pt_on_curve, junk, &tangent1 );
      if( closest_coedge->sense() == CUBIT_REVERSED )
        tangent1 *= -1.0;
        
      inside = (inside_of_curve(tangent1, pt_on_curve, pt_on_surf, normal) == CUBIT_TRUE);
    }
  }
    // **** Special case: closest point on loop at vertex ****
  else
  {
    CURVE* curve1_ptr = closest_coedge->get_curve();
    CURVE* curve2_ptr =   other_coedge->get_curve();
    curve1_ptr->closest_point( pt_on_curve, junk, &tangent1 );
    curve2_ptr->closest_point( pt_on_curve, junk, &tangent2 );

    if( closest_coedge->sense() == CUBIT_REVERSED )
      tangent1 *= -1.0;
    if(   other_coedge->sense() == CUBIT_REVERSED )
      tangent2 *= -1.0;

    surface->closest_point( pt_on_curve, NULL, &coe_normal );
    
    coe_cross = tangent1 * tangent2;
    double sum  = (coe_cross + coe_normal).length_squared();
    double diff = (coe_cross - coe_normal).length_squared();
    
    CubitBoolean inside1 = 
      inside_of_curve( tangent1, pt_on_curve, pt_on_surf, normal );
    CubitBoolean inside2 = 
      inside_of_curve( tangent2, pt_on_curve, pt_on_surf, normal );
        
    if( (sum > diff) || (diff < CUBIT_DBL_MIN) ) 
      //the common vertex is at a convexity
    {
      //the point must be inside of both coedges
      inside = inside1 && inside2;
    }
    else //the common vertex is at a concavity
    {
      //the point must be inside of one of the coedges
      inside = inside1 || inside2;
    }
  }
        
  result_pt = inside ? pt_on_surf : pt_on_curve;
  if( containment )
  {
    if( (pt_on_surf-pt_on_curve).length_squared() < (GEOMETRY_RESABS*GEOMETRY_RESABS) )
      *containment = CUBIT_PNT_BOUNDARY;
    else 
      *containment = inside ? CUBIT_PNT_INSIDE : CUBIT_PNT_OUTSIDE;
  }

  return CUBIT_SUCCESS;
}  

template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
CubitStatus VGLoopTool<SURFACE,LOOP,COEDGE,CURVE,POINT>::
closest_surface_coedge( SURFACE* surface, const CubitVector& from_pt,
                        COEDGE*& coedge1, COEDGE*& coedge2,
                        CubitVector* pt_on_curve )
{
  coedge1 = coedge2 = 0;
  LOOP* loop = 0;
  COEDGE* coedge;
  const double ressqr = CUBIT_RESABS * CUBIT_RESABS;


  
  //Find the first bounding box that the point is within, or if the 
  // point isn't in any bounding box, the closest one
  closest_dst = CUBIT_DBL_MAX;
  while( loop = surface->next_loop( loop ) )
  {
    coedge = loop->first_coedge();
    do
    {
      if( odd_coedge_count( coedge ) )
      {
        double dist_sqr = 
          coedge->get_curve()->bounding_box().distance_squared(from_pt);
        if( dist_sqr < closest_dst )
        {
          coedge1 = coedge;
          closest_dst = dist_sqr;
          if( dist_sqr < ressqr )
            break;
        }
      }
      
      coedge = coedge->next();
    } while( coedge != loop->first_coedge() );
    
    if( closest_dst < ressqr )
      break;
  }

  if( !coedge1 )
    return CUBIT_FAILURE;
  
  //Start by assuming that the curve we found above is the closest
  CubitVector closest_pt;
  coedge1->get_curve()->closest_point_trimmed( from_pt, closest_pt );
  closest_dst = ( from_pt - closest_pt ).length_squared();
  
  //Now look for a closer curve
  while( loop = surface->next_loop( loop ) )
  {
    coedge = loop->first_coedge();
    do
    {
      if( odd_coedge_count( coedge ) )
      {
        if( coedge->get_curve()->bounding_box().distance_squared(from_pt) <=
          closest_dst+GEOMETRY_RESABS )
        {
          CubitVector current_pt;
          coedge->get_curve()->closest_point_trimmed( from_pt, current_pt );
          double current_dst = (from_pt - current_pt).length_squared();

          if( (closest_pt = current_pt).length_squared() < ressqr )
          {
            if( current_dst < closest_dst )
            {
              closest_dst = current_dst;
              closest_pt = current_pt;
              coedge2 = coedge;
            }
            else
            {
              coedge2 = coedge;
            }
          }
          else if( current_dst < closest_dst )
          {
            coedge1 = coedge;
            closest_dst = current_dst;
            closest_pt = current_pt;
            coedge2 = 0;
          }
        }
      }
      coedge = coedge->next();
    } while( coedge != loop->first_coedge() );
  }
  
  
  //make sure we have things in the correct order
  if( coedge2 )
  {
    POINT* common = coedge1->get_curve()->common( coedge2->get_curve() );
    if( ! common ) 
      coedge2 = 0;
    else if( coedge1->start_point() == coedge2->end_point() )
    {
      COEDGE* tmp = coedge1;
      coedge1 = coedge2;
      coedge2 = tmp;
    }
  }

  if( pt_on_curve ) *pt_on_curve = closest_pt;
  return CUBIT_SUCCESS;
}  
*/
