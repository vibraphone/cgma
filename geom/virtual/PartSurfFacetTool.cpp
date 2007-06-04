#include "PartSurfFacetTool.hpp"
#include "CubitMessage.hpp"
#include "DLIList.hpp"
#include "GeometryQueryEngine.hpp"

#include "PartitionSurface.hpp"
#include "PartitionCurve.hpp"
#include "PartitionPoint.hpp"
#include "TDVGFacetOwner.hpp"

#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"

#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "GfxDebug.hpp"
#include "GMem.hpp"

#undef PART_SURF_REFACET

void PartSurfFacetTool::validate_facets( PartitionSurface* mySurface )
{
  int i, j, k;
  DLIList<CubitFacetData*> facets;
  DLIList<CubitFacet*> pt_facets, edge_facets;
  DLIList<CubitFacetEdge*> pt_edges;
  mySurface->get_facet_data(facets);
  for( i = facets.size(); i--; )
  {
    CubitFacetData* facet = facets.get_and_step();
    for ( j = 0; j < 3; j++ )
    {
      CubitFacetEdge* edge = facet->edge(j);
      if ( !edge )
      {
        PRINT_ERROR("Facet with NULL edge.\n");
        continue;
      }
      
      if ( PartitionEntity* ent = TDVGFacetOwner::get(edge) )
      {
        PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent);
        if (!curve)
        {
          PRINT_ERROR("Facet edge owned by non-curve %s %p\n",
            typeid(*ent).name(), (void*)ent);
          GfxDebug::draw_line( edge->point(0)->coordinates(),
                               edge->point(1)->coordinates(),
                               CUBIT_RED );
          continue;
        }
        
        int edge_use_count = 0;
        for ( k = 0; k < edge->num_adj_facets(); k++ )
          if ( (const PartitionEntity*)TDVGFacetOwner::get(edge->adj_facet(k)) == mySurface )
            edge_use_count++;
        assert(edge_use_count);
        
        int curve_use_count = 0;
        int coedge_count = 0;
        PartitionCoEdge* coedge = 0;
        while (( coedge = curve->next_coedge(coedge) ))
        {
          if (coedge->get_loop())
          {
            coedge_count++;
            if (coedge->get_loop()->get_surface() == mySurface)
              curve_use_count++;
          }
        }
          
        RefEdge* curve_owner = dynamic_cast<RefEdge*>(curve->topology_entity());
        if ( coedge_count != edge->num_adj_facets() )
        {
          PRINT_ERROR("Curve %p (RefEdge %d): Curve in %d partition surfaces "
                      "has facet edge in %d facets.\n", (void*)curve,
                      curve_owner ? curve_owner->id() : 0,
                      coedge_count, edge->num_adj_facets() );
          GfxDebug::draw_line( edge->point(0)->coordinates(),
                               edge->point(1)->coordinates(),
                               CUBIT_RED );
        }
        
        if ( curve_use_count != edge_use_count )
        {
          PRINT_ERROR("Curve %p (RefEdge %d): Curve used in %d coedges of "
                      "this surface exists in %d facets of this surface.\n", 
                      (void*)curve,
                      curve_owner ? curve_owner->id() : 0,
                      curve_use_count, edge_use_count );
           GfxDebug::draw_line( edge->point(0)->coordinates(),
                                edge->point(1)->coordinates(),
                                CUBIT_RED );
        }
      }
      else  // FacetEdge has no owner
      {
        if (edge->num_adj_facets() != 2)
        {
          PRINT_ERROR("Interior facet edge is %d-valent.\n", edge->num_adj_facets());
          GfxDebug::draw_line( edge->point(0)->coordinates(),
                               edge->point(1)->coordinates(),
                               CUBIT_RED );
        }
      }
      
      bool found = false;
      for ( k = 0; k < edge->num_adj_facets(); k++ )
      {
        CubitFacet* edge_facet = edge->adj_facet(k);
        if( edge_facet == facet )
          found = true;
        
        int index = edge_facet->edge_index(edge);
        if( index < 0 )
        {
          PRINT_ERROR("edge->facet link missing facet->edge link.\n");
          continue;
        }
        CubitPoint* fp1 = edge_facet->point((index+1)%3);
        CubitPoint* fp2 = edge_facet->point((index+2)%3);
        int sense = 1;
        if( fp1 == edge->point(1) && fp2 == edge->point(0) )
          sense = -1;
        else if( fp1 != edge->point(0) || fp2 != edge->point(1) )
        {
          PRINT_ERROR("Facet adjacent to edge does not contain edge end points.\n");
          GfxDebug::draw_line( edge->point(0)->coordinates(),
                               edge->point(1)->coordinates(),
                               CUBIT_RED );
          continue;
        }
        if ( sense != edge_facet->edge_use(index) )
        {
          //PRINT_ERROR("Facet has incorrect edge use\n");
        }
      }
      
      if( !found )
        PRINT_ERROR("facet->edge link missing edge->facet link.\n");
    }
    
    for ( j = 0; j < 3; j++ )
    {
      CubitPoint* point = facet->point(j);
      pt_facets.clean_out();
      point->facets(pt_facets);
      bool found = false;
      for ( k = pt_facets.size(); k--; )
      {
        CubitFacet* pt_facet = pt_facets.get_and_step();
        if ( pt_facet == facet )
          found = true;
        else if( pt_facet->point_index(point) < 0 )
          PRINT_ERROR("point->facet link missing facet->point link.\n");
      }
      if( ! found )
        PRINT_ERROR("facet->point link missing point->facet link.\n");
        
      point->edges(pt_edges);
      while (pt_edges.size())
      {
        CubitFacetEdge* edge = pt_edges.pop();
        edge_facets.clean_out();
        edge->facets(edge_facets);
        int count = 0;
        while (edge_facets.size())
          if (TDVGFacetOwner::get(edge_facets.pop()) == mySurface)
            count++;
        
        bool draw = false;
        if (count == 0) // not in this face
         ;
        else if (PartitionEntity* ent = TDVGFacetOwner::get(edge))
        {
          PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent);
          if (curve->is_nonmanifold(mySurface))
          {
            if (count != 2)
            {
              draw = true;
              PRINT_ERROR("Edge on non-manifold curve contained in %d facets.\n", count);
            }
          } 
          else if (count != 1)
          {
            draw = true; 
            PRINT_ERROR("Boundary edge contained in %d facets.\n", count); 
          }
        }
        else if (count != 2)
        { 
          draw = true; 
          PRINT_ERROR("Interior edge contained in %d facets.\n", count); 
        }
      
        if (draw) {
          CubitVector start = edge->point(0)->coordinates();
          CubitVector step = 0.05 * (edge->point(1)->coordinates() - start);
          for (int k = 0; k < 20; k++)
            GfxDebug::draw_point( start + k * step, CUBIT_RED );
          GfxDebug::flush();
        }
      }
    }
  }
  
  
    // Check loops for contiguous chains of facet edges.
  DLIList<CubitFacetEdgeData*> curve_edges, loop_edges;
  PartitionLoop* loop = 0;
  while ((loop = mySurface->next_loop(loop)))
  {
    loop_edges.clean_out();
    PartitionCoEdge* coedge = loop->first_coedge();
    do {
      PartitionCurve* curve = coedge->get_curve();
      if (!curve->is_nonmanifold(mySurface))
      {
        curve_edges.clean_out();
        curve->get_facet_data(curve_edges);
        if (coedge->sense() == CUBIT_REVERSED)
          curve_edges.reverse();
        loop_edges += curve_edges;
      }
      coedge = loop->next_coedge(coedge);
    } while (coedge != loop->first_coedge());
  
    loop_edges.last();
    CubitFacetEdgeData* prev = loop_edges.get_and_step();
    PartitionCurve* prev_curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(prev));
    for (int i = 0; i < loop_edges.size(); i++)
    {
      CubitFacetEdgeData* edge = loop_edges.get_and_step();
      CubitPoint* shared = edge->shared_point(prev);
      PartitionCurve* curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edge));
      
      if (!prev_curve || !curve)
      {
        PRINT_ERROR("Missing facet owner on boundary curve.\n");
        prev_curve = curve;
        prev = edge;
        continue;
      }
      
      if (!shared)
      {
        PRINT_ERROR("Broken facet-edge chain at curve %p (RefEdge %d):\n"
                    "      start of segment (%f,%f,%f)->(%f,%f,%f)\n",
                    curve, curve && dynamic_cast<RefEdge*>(curve->topology_entity())?
                    dynamic_cast<RefEdge*>(curve->topology_entity())->id() : 0,
                    edge->point(0)->coordinates().x(), 
                    edge->point(0)->coordinates().y(),
                    edge->point(0)->coordinates().z(),
                    edge->point(1)->coordinates().x(), 
                    edge->point(1)->coordinates().y(),
                    edge->point(1)->coordinates().z());
      }
      else if(PartitionPoint *pt = dynamic_cast<PartitionPoint*>(TDVGFacetOwner::get(shared)))
      {
        if( (prev_curve == curve) &&
             !(pt == curve->start_point() && pt == curve->end_point()) )
        {
          PRINT_ERROR("Edges at point %p (vertex %d) (%f,%f,%f) owned\n"
                      "      by the same curve: %p (RefEdge %d)\n",
            pt, dynamic_cast<RefVertex*>(pt->topology_entity()) ?
            dynamic_cast<RefVertex*>(pt->topology_entity())->id() : 0,
            pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z(),
            curve, curve && dynamic_cast<RefEdge*>(curve->topology_entity())?
            dynamic_cast<RefEdge*>(curve->topology_entity())->id() : 0);
        } 
        else if (!curve->other_point(pt))
        {
          PRINT_ERROR("Edge owner curve %p (RefEdge %d) at point %p \n"
                      "      (vertex %d) (%f,%f,%f) doesn't match topology\n",
            curve, curve && dynamic_cast<RefEdge*>(curve->topology_entity())?
            dynamic_cast<RefEdge*>(curve->topology_entity())->id() : 0,
            pt, dynamic_cast<RefVertex*>(pt->topology_entity()) ?
            dynamic_cast<RefVertex*>(pt->topology_entity())->id() : 0,
            pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z());
        }
        else if(!prev_curve->other_point(pt))
        {
          PRINT_ERROR("Edge owner curve %p (RefEdge %d) at point %p \n"
                      "      (vertex %d) (%f,%f,%f) doesn't match topology\n",
            prev_curve, prev_curve && dynamic_cast<RefEdge*>(prev_curve->topology_entity())?
            dynamic_cast<RefEdge*>(prev_curve->topology_entity())->id() : 0,
            pt, dynamic_cast<RefVertex*>(pt->topology_entity()) ?
            dynamic_cast<RefVertex*>(pt->topology_entity())->id() : 0,
            pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z());
        }
      }
      else if (prev_curve != curve)
      {
        PRINT_ERROR("Missing point-owner at  transition between\n"
                  "      curve %p (RefEdge %d) and curve %p (RefEdge %d):\n"
                  "      at location (%f,%f,%f)\n",
                  curve, curve && dynamic_cast<RefEdge*>(curve->topology_entity())?
                  dynamic_cast<RefEdge*>(curve->topology_entity())->id() : 0,
                  prev_curve, prev_curve && dynamic_cast<RefEdge*>(prev_curve->topology_entity())?
                  dynamic_cast<RefEdge*>(prev_curve->topology_entity())->id() : 0,
                  shared->coordinates().x(), 
                  shared->coordinates().y(),
                  shared->coordinates().z() );
      } 
      
      prev = edge;
      prev_curve = curve;
    }
  }
}
   

//-------------------------------------------------------------------------
// Purpose       : Find closest point on passed facet
//
// Special Notes : static
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/28/03
//-------------------------------------------------------------------------
void PartSurfFacetTool::closest_pt_on_facet( CubitFacet* facet,
                                             const CubitVector& p,
                                             CubitVector& result )
{
    // Get triangle vertices
  CubitVector p0(facet->point(0)->coordinates());
  CubitVector p1(facet->point(1)->coordinates());
  CubitVector p2(facet->point(2)->coordinates());

/*
  Algorithm from:
    "Distance Between Point and Triangle in 3D"
    David Eberly
    Magic Software, Inc.
    Sept. 28, 1999

  Use barycentric coordinates.  Coordinates are
  calculated in the range [0,det] rather than [0,1]
  to avoid the fp division entirely where it can
  be avoided.

    ^v*t                                 
 \R2|                 
  \ |                 
   \|                 
    *p2               
    |\                
    | \     R1        
 R3 |  \              
    |   \             
    | R0 \            
    |     \p1         
 ---*------*--->u*s   
    |p0     \  R6     
 R4 |   R5   \        
    |         \u+v=det  
*/


  CubitVector s(p1 - p0);  // the u (or s) axis in parameterized space
  CubitVector t(p2 - p0);  // the v (or t) axis in parameterized space
  CubitVector d(p0 - p);   // vector from input position to corner at (u,v) = (0,0)
    // Pre-calculate all the dot products we need
    // Name the dot product of vectors 'a' and 'b' as 'ab'
  double ss = s.length_squared();
  double st = s % t;
  double tt = t.length_squared();
  double sd = s % d;
  double td = t % d;
    // Calculate barycentric coordinates in the range [0,det]
  double det = ss*tt - st*st;
  double u = st*td - tt*sd;
  double v = st*sd - ss*td;
  
    // Big tree of conditionals to determine which of the 
    // regions in the above diagram the projection of
    // the point into the plane lies in.
  if ( u+v < det )
  {
    if ( u < 0 )
    {
      if( v < 0 )    // Region 4 
      {
        if ( sd < 0 )
        {
          if ( -sd > ss )
            result = p1;
          else
            result = p0 + (-sd/ss) * s;
        }
        else if ( td < 0 )
        {
          if ( -td > tt )
            result = p2;
          else
            result = p0 + (-td/tt) * t;
        }
        else
        {
          result = p0;
        }
      }
      else           // Region 3 (Edge p2-p0, u->0)
      {
        if ( td > 0 )
          result = p0;
        else if ( -td > tt )
          result = p2;
        else
          result = p0 + (-td/tt) * t;
      }
    }
    else if ( v < 0) // Region 5 (Edge p0-p1, v->0)
    {
      if ( sd > 0 )
        result = p0;
      else if ( -sd > ss )
        result = p1;
      else 
        result = p0 + (-sd/ss) * s;
    }
    else             // Region 0 (Interior)
    {
      result = p0 + (1.0/det) * (u*s + v*t);
    }
  }
  else if ( u < 0 )  // Region 2 
  {
    double num = tt + td - st - sd;
    if ( num > 0.0 )       
    {
      double den = ss - 2*st + tt;
      if ( num >= den )    // (Point p1)
        result = p1;
      else                 // (Edge p1-p2)
      {
        u = num / den;
        result = p0 + u*s + (1-u)*t;
      }
    } 
    else if ( td >= 0 )    // (Point p0)
      result = p0;
    else if ( tt+td <= 0 ) // (Point p2)
      result = p2;
    else                   // (Edge p2-p0)
      result = p0 + (-td/tt)*t;
  }
  else if ( v < 0 )  // Region 6 
  {
    double num = tt + td - st - sd;
    double den = ss - 2*st + tt;
    if (num < den)       
    {
      if ( num < 0 )    // (Point p1)
        result = p2;
      else                 // (Edge p1-p2)
      {
        u = num / den;
        result = p0 + u*s + (1-u)*t;
      }
    } 
    else if ( sd >= 0 )    // (Point p0)
      result = p0;
    else if ( ss+sd <= 0 ) // (Point p1)
      result = p1;
    else                   // (Edge p0-p1)
      result = p0 + (-sd/ss)*s;
  }
  else               // Region 1 (Edge p1-p2, u+v->1)
  {
    double num = tt + td - st - sd;
    if ( num <= 0 )
      result = p2;
    else
    {
      double den = ss - 2*st + tt;
      if ( num >= den )
        result = p1;
      else
      {
        u = num/den;
        result = p0 + u*s + (1-u)*t;
      }
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Local method for debug output
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/03
//-------------------------------------------------------------------------
static void draw_edges( DLIList<CubitFacetEdge*>& edges, int edge_color = 0, 
                 bool label_edges = false, bool draw_points = false,
                 int point_color = 0, bool label_points = false )
{
  char buffer[2*sizeof(void*)+3];

  if (edge_color == 0)
    edge_color = CUBIT_WHITE;
  if (point_color == 0)
    point_color = CUBIT_WHITE;
  
  for (int i = edges.size(); i--; )
  {
    CubitFacetEdge* edge = edges.get_and_step();
    CubitVector start = edge->point(0)->coordinates();
    CubitVector end = edge->point(1)->coordinates();
    GfxDebug::draw_line(start, end, edge_color );
    
    if (label_edges)
    {
      CubitVector mid = 0.5 * (start + end);
      sprintf(buffer, "%p", edge);
      float x = (float)mid.x();
      float y = (float)mid.y();
      float z = (float)mid.z();
      GfxDebug::draw_label( buffer, x, y, z, edge_color );
    }
      
    for ( int j = 0; j < 2; j++ )
    {
      CubitPoint* point = edge->point(j);
      float x = (float)point->coordinates().x();
      float y = (float)point->coordinates().y();
      float z = (float)point->coordinates().z();
      
      if (draw_points)
      {
        GfxDebug::draw_point( x, y, z, point_color );
      }
      if (label_points)
      {
        sprintf(buffer, "%p", point);
        GfxDebug::draw_label( buffer, x, y, z, point_color );
      }
    }
  }
  GfxDebug::flush();
}

    
CubitStatus PartSurfFacetTool::init_facet_data( 
                                 DLIList<CubitFacetData*>& facets )
{
  assert(!mySurface->has_facets());
  
  int i;
  CubitStatus rval;
  DLIList<CubitPoint*> boundary_points, interior_points;
  DLIList<CubitFacetEdge*> boundary_edges, interior_edges;
  DLIList<PartitionPoint*> geom_points, interior_geom_points;
  
    // Get facet data
  rval = get_facet_points_and_edges( facets, boundary_points, 
    interior_points, boundary_edges, interior_edges );
  if (!rval)
  {
    PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__ );
    return CUBIT_FAILURE;
  }

    // Get real points on surface boundary
  mySurface->get_points(geom_points);
  geom_points.last();
  for (i = geom_points.size(); i--; )
    if (!geom_points.step_and_get()->real_point())
      geom_points.change_to(0);
  geom_points.remove_all_with_value(0);
  
    // Group geometric points into boundary and interior sets
  for (i = geom_points.size(); i--; )
  {
    PartitionPoint* point = geom_points.step_and_get();
    bool boundary = false;
    PartitionCurve* curve = 0;
    while ( (curve = point->next_curve(curve)) )
      if (curve->is_in_surface(mySurface,true))
        boundary = true;
    if (!boundary)
    {
      geom_points.change_to(0);
      interior_geom_points.append(point);
    }
  }
  geom_points.remove_all_with_value(0);
  
    // Associate each boundary geometry point with the closest
    // boundary facet point.
  if(boundary_points.size() == 0 && geom_points.size() != 0)
    rval = CUBIT_FAILURE;
  else
    rval = associate_points( boundary_points, geom_points );

  if (!rval)
    return CUBIT_FAILURE;
 
    // Associate each interior geometry point with the closest 
    // interior facet point.
  rval = associate_points( interior_points, interior_geom_points );
  if (!rval)
    return CUBIT_FAILURE;

  if (DEBUG_FLAG(145))
  {
    draw_edges( boundary_edges, CUBIT_ORANGE, true, false, CUBIT_WHITE, true );
    draw_edges( interior_edges, CUBIT_WHITE );
  }
  
  
    // Populate sets
  boundary_set.clear();
  interior_set.clear();
  for (i = boundary_edges.size(); i--; )
    boundary_set.insert(boundary_edges.step_and_get());
  for (i = interior_edges.size(); i--; )
    interior_set.insert(interior_edges.step_and_get());
 
    // For each loop on the surface...
  PartitionLoop* loop = 0;
  DLIList<PartitionCurve*> curve_set;
  DLIList<CubitFacetEdge*> point_edges, edge_set_1, edge_set_2;
  DLIList<PartitionCurve*> nonmanifold_stack;
  while ( (loop = mySurface->next_loop(loop)) )
  {
      // Find coedge beginning at real vertex (must be at least one)
    PartitionCoEdge* coedge = loop->first_coedge();
    while (!coedge->start_point()->real_point())
    { 
      coedge = coedge->next();
      assert(coedge != loop->first_coedge());
    }
    
    PartitionCoEdge* first_coedge = coedge; // so we know when to stop
    
     
    do // Loop until all curves in loop have been handled (back to first_coedge)
    {
        // Keep track of start and end
      PartitionCurve* first_curve = coedge->get_curve();
      PartitionPoint* start_point = coedge->start_point();
      PartitionPoint* end_point = coedge->end_point();

        // Get list of curves until next real vertex
      curve_set.clean_out();
      curve_set.append(first_curve);
      while (!coedge->end_point()->real_point())
      {
        coedge = coedge->next();
        curve_set.append(coedge->get_curve());
        end_point = coedge->end_point();
        assert(&first_curve->sub_entity_set() == &coedge->get_curve()->sub_entity_set());
      }
      
      if (DEBUG_FLAG(145))
      {
        curve_set.reset();
        for (i = curve_set.size(); i--; )
        {
          GMem gmem;
          int junk;
          PartitionCurve* c = curve_set.get_and_step();
          c->get_geometry_query_engine()->get_graphics( c, junk, &gmem );
          GfxDebug::draw_polyline(gmem.point_list(), gmem.pointListCount, CUBIT_RED );
        }
        GfxDebug::flush();
      }

        // Special case for non-manifold curves
      if (first_curve->is_nonmanifold(mySurface))
      {
          // Don't do the same non-manifold curve twice
        curve_set.reset();
        nonmanifold_stack.last();
        if (nonmanifold_stack.size() &&
            curve_set.get() == nonmanifold_stack.get())
        {
          for (i = curve_set.size(); i--; )
          {
            PartitionCurve* curv1 = curve_set.get_and_step();
            PartitionCurve* curv2 = nonmanifold_stack.pop();
            assert(curv1 == curv2);
          }
        }
        else
        {
          for (i = curve_set.size(); i--; )
            nonmanifold_stack.append(curve_set.get_and_step());
        
          if (!seam_nonmanifold_curves(curve_set, facets))
          {
            return CUBIT_FAILURE;
          }
        }
        coedge = coedge->next();
        if (coedge != first_coedge)
          continue;
        else
          break;
      }

        // Get facet edges adjacent to first point
      point_edges.clean_out();
      start_point->facet_point()->edges(point_edges);

        // For each boundary edge on the point, look for a chain
        // of connected boundary edges.
      edge_set_1.clean_out();
      edge_set_2.clean_out();
      while(point_edges.size())
      {
        CubitFacetEdge* edge = point_edges.pop();
        if (boundary_set.count(edge) == 0)  // skip non-boundary edges
          continue;

        edge_set_2.clean_out();
        if (!get_boundary_chain(start_point->facet_point(), edge, 
                                  end_point->facet_point(), edge_set_2))
          continue; // skip edges that don't begin a chain that reaches the end vertex

          // If only one chain found so far, store the chain and 
          // continue to the next edge around the point
        if (!edge_set_1.size())
        {
          edge_set_1 = edge_set_2;
          edge_set_2.clean_out();
          continue;
        }

          // If we reached this point, handle special case where
          // loop is composed of two real curves and we need a 
          // geometric check to determine which half of the facet edge
          // loop goes with which curve.
        edge_set_1.reset();
        edge_set_2.reset();
        CubitVector set_1_mid = edge_set_1.get()->center();
        CubitVector set_2_mid = edge_set_2.get()->center();
        CubitVector set_1_closest, set_2_closest;
        first_curve->closest_point_trimmed(set_1_mid, set_1_closest);
        first_curve->closest_point_trimmed(set_2_mid, set_2_closest);

        set_1_closest -= set_1_mid;
        set_2_closest -= set_2_mid;
        if (set_2_closest.length_squared() < set_1_closest.length_squared())
          edge_set_1 = edge_set_2;
      } // while(point_edges.size())

        // We have the chain of edges to associate with the curve
        // set.  Seam with adjacent surface facettings (if any)
      if (DEBUG_FLAG(145))
        draw_edges( edge_set_1, CUBIT_RED, true );
      
      if (!seam_curves( curve_set, edge_set_1, facets ))
        return CUBIT_FAILURE;
      
      coedge = coedge->next();
    } while (coedge != first_coedge);
    
  } // while (loop)
  
    // attach facets to surface
  mySurface->set_facet_data(facets);
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Helper function for curve seaming functions.
//                 Check that input curve list are partitions of the
//                 same real curve and are all the partitions of that curve.
//                 Returns the real curve and the PartitionPoints that
//                 are the start and end of the real curve.
//
// Special Notes : Assumes the passed curves are either in order or in
//                 the reverse order.  If the reverse order, the passed
//                 list will be reversed.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/08/03
//-------------------------------------------------------------------------
Curve* PartSurfFacetTool::get_real_curve( DLIList<PartitionCurve*>& curve_list,
                                          PartitionPoint*& start_point,
                                          PartitionPoint*& end_point )
{
  SubEntitySet* set = &(curve_list.get()->sub_entity_set());
  DLIList<PartitionEntity*> tmp_set(curve_list.size());
  set->get_sub_entities(tmp_set);
  if (tmp_set.size() != curve_list.size())
    return 0;
  
  tmp_set.reset();
  curve_list.reset();
  if (tmp_set.get() != curve_list.get())
  {
    curve_list.reverse();
    curve_list.reset();
  }
  
  for (int i = curve_list.size(); i--; )
    if (curve_list.get_and_step() != tmp_set.get_and_step())
      return 0;
   
    // Get real curve
  Curve* real_curve = dynamic_cast<Curve*>(set->get_entity());
  assert(!!real_curve);
  
  DLIList<TopologyBridge*> point_bridges(2);
  real_curve->get_children(point_bridges, false, set->get_owner_layer());
  point_bridges.reset();
  start_point = dynamic_cast<PartitionPoint*>(point_bridges.get());
  end_point = dynamic_cast<PartitionPoint*>(point_bridges.next());
  assert(start_point && end_point);
  
  return real_curve;
}

CubitStatus PartSurfFacetTool::get_boundary_chain( 
                                      CubitPoint* start_point,
                                      CubitFacetEdge* start_edge,
                                      CubitPoint* end_point,
                                      DLIList<CubitFacetEdge*>& result_list )
{
  DLIList<CubitFacetEdge*> point_edges;
  
  //assert(start_edge->marked());
  assert(result_list.size() == 0);
  CubitPoint* point = start_point;
  CubitFacetEdge* edge = start_edge;
  while (true)
  {
    result_list.append(edge);
    point = edge->other_point(point);
    if (TDVGFacetOwner::get(point))
      return point == end_point ? CUBIT_SUCCESS : CUBIT_FAILURE;
    
    point->edges(point_edges);
    for (int i = point_edges.size(); i--; )
    {
      CubitFacetEdge* point_edge = point_edges.step_and_get();
      //if (point_edge->marked() != 1 || point_edge == edge)
      if (boundary_set.count(point_edge) == 0 || point_edge == edge)
        point_edges.change_to(0);
    }
    point_edges.remove_all_with_value(0);
    
    if (point_edges.size() != 1)
      return CUBIT_FAILURE;
    
    edge = point_edges.remove();
  }
  
    // unreachable
  return CUBIT_SUCCESS;
}
      
//-------------------------------------------------------------------------
// Purpose       : Associate facet edges in the interior of a surface
//                 facetting with the set of curve partitions of some
//                 real, non-manifold curve in the surface.
//
// Special Notes : Assumes edges internal to the surface facetting
//                 have been marked with a '2' by the caller.  Clears
//                 marks on consumed edges.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/08/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::seam_nonmanifold_curves(
                                        DLIList<PartitionCurve*>& curve_list,
                                        DLIList<CubitFacetData*>& facet_list )
{
    // All the curves in the attached coedge_list should belong
    // to the same real curve and thus have a single continuous 
    // parameterization over the entire set of curves.  Verify
    // this now.
  PartitionPoint *start_point, *end_point;
  Curve* real_curve = get_real_curve(curve_list, start_point, end_point);
  if (!real_curve)
  {
    PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__ );
    assert(false);
    return CUBIT_FAILURE;
  }
  
    // Assume caller has marked all interior edges with a '2'.
    // Clear marks as we go to ensure we don't loop forever.
    // Find chain of curves closest to real curve and ending
    // at the end facet point.
  CubitPoint* const start_vtx = start_point->facet_point();
  CubitPoint* const   end_vtx =   end_point->facet_point();
  assert(start_vtx&&end_vtx);
  DLIList<CubitFacetEdge*> vtx_edges, seam_list;
  CubitPoint* vtx = start_vtx;
  do
  {
    vtx->edges(vtx_edges);
    CubitFacetEdge* edge = 0;
    double shortest_dist_sqr = CUBIT_DBL_MAX;
    while( vtx_edges.size() )
    {
      CubitFacetEdge* curr = vtx_edges.pop();
      //if (curr->marked() != 2 )
      if (interior_set.count(curr) == 0)
        continue;
      
      if (curr->other_point(vtx) == end_vtx)
      { 
        edge = curr;
        break;
      }
      
      CubitVector closest, tangent;
      CubitVector start(vtx->coordinates());
      CubitVector end(curr->other_point(vtx)->coordinates());
      real_curve->closest_point(0.5*(start+end), closest, &tangent);
      if (tangent % (end-start) < 0.0)
        continue;
      
      real_curve->closest_point(end, closest);
      closest -= end;
      double dist_sqr = closest.length_squared();
      if (dist_sqr < shortest_dist_sqr)
      {
        edge = curr;
        shortest_dist_sqr = dist_sqr;
      }
    }
    
    if (!edge)
    {
      PRINT_ERROR("Internal error at %s:%d:\n"
                  "Failed to associate facet edges with non-manofold curves.\n",
                  __FILE__, __LINE__ );
      return CUBIT_FAILURE;
    }
    if (DEBUG_FLAG(145))
    {
      GfxDebug::draw_facet_edge( edge, CUBIT_MAGENTA );
      GfxDebug::flush();
    }
    
    seam_list.append(edge);
    //edge->marked(0);
    vtx = edge->other_point(vtx);
  } while (vtx != end_vtx);
  
  return seam_curves( curve_list, seam_list, facet_list );
}


CubitStatus PartSurfFacetTool::split_edge( CubitFacetEdge* old_edge,
                                          const CubitVector& position,
                                          CubitFacet* edge_facet,
                                          CubitPoint*& new_point,
                                          CubitFacetEdge*& new_edge,
                                          CubitFacet*& new_facet )
{
  CubitVector v1, v2;
  int i = 0, junk;
  new_point = 0;
  new_edge = 0;
  new_facet = 0;

  CubitFacet* split_facet = edge_facet;
  if( !split_facet )
    split_facet = old_edge->adj_facet(0);
  else if( edge_facet->edge_index(old_edge) < 0 )
    { assert(0); return CUBIT_FAILURE; }

  CubitPoint* pt1 = old_edge->point(0);
  CubitPoint* pt2 = old_edge->point(1);

  v1 = position - pt1->coordinates();
  v2 = pt2->coordinates() - position;
  assert( v1.length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS );
  assert( v2.length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS );
  double dot = v1 % v2;
  assert( dot > 0 ); // projection of new point onto line lies outside edge.
//  double cos_sqr = (dot * dot) / (v1.length_squared() * v2.length_squared());
//  assert( cos_sqr > 0.99240387650610407 ); // angle less than 5 degrees


  CubitPoint* new_pt = split_facet->split_edge( pt1, pt2, position );
  new_edge = new_pt->shared_edge(pt2);
  if( !new_edge ) { assert(0); return CUBIT_FAILURE; }

    // find new facet and update for other split facets
  CubitFacet *facet, *other_facet;
  PartitionEntity* facet_owner;
  while ( (facet = old_edge->adj_facet(i++)) ) {
    CubitPoint* pt3 = facet->point( facet->edge_index(pt1,new_pt,junk) );
    if( !pt3->shared_edge( new_pt ) )
      new CubitFacetEdgeData( pt3, new_pt );
      
    other_facet = facet->shared_facet( new_pt, pt3 );
    if( !other_facet ) { assert(0); continue; }
    if ( facet == edge_facet )
      new_facet = other_facet;
    else if( (facet_owner = TDVGFacetOwner::get(facet)) )
      facet_owner->notify_split( facet, other_facet );
  }
  if( (facet_owner = TDVGFacetOwner::get( old_edge )) )
      facet_owner->notify_split( old_edge, new_edge );

  if( edge_facet )
    assert(!!new_facet);
  else 
    new_facet = 0;
    
  new_point = new_pt;
  return CUBIT_SUCCESS;
}

CubitStatus PartSurfFacetTool::collapse_edge( CubitPoint* keep,
                                             CubitPoint* dead,
                                             DLIList<CubitFacetData*>* unowned )
{
  int i;
  
  CubitPointData* keep_pt = dynamic_cast<CubitPointData*>(keep);
  CubitPointData* dead_pt = dynamic_cast<CubitPointData*>(dead);
  CubitFacetEdge* dead_edge;
  assert(keep_pt&&dead_pt);
  
    // Cannot proceed if dead_pt is owned by a PartitionPoint.
  if ( TDVGFacetOwner::get(dead_pt) )
    return CUBIT_FAILURE;
  
    // Get list of facets to be destroyed when edge is collapsed
  DLIList<CubitFacet*> dead_facets;
  keep->shared_facets(dead_pt,dead_facets);
  DLIList<PartitionSurface*> dead_facet_owners(dead_facets.size());
  DLIList<CubitFacetData*> dead_facet_ds(dead_facets.size());
  CAST_LIST(dead_facets, dead_facet_ds, CubitFacetData);
  
    // Find other edges to be destroyed when edge is collapsed.
    // cannot proceed if any of them belong to a PartitionCurve.
    // Also, find owners of each facet to be destroyed.
  dead_facets.reset();
  for ( i = dead_facets.size(); i--; )
  {
    CubitFacet* facet = dead_facets.get_and_step();
    int keep_index = facet->point_index(keep);
    int dead_index = facet->point_index(dead);
    int othr_index = (keep_index + 1) % 3;
    if ( othr_index == dead_index )
      othr_index = (keep_index + 2) % 3;
    assert( keep_index >= 0 && dead_index >= 0 );
    dead_edge = dead_pt->shared_edge(facet->point(othr_index));
    if ( TDVGFacetOwner::get(dead_edge) )
      return CUBIT_FAILURE;
      
    PartitionSurface* surf = dynamic_cast<PartitionSurface*>(TDVGFacetOwner::get(facet));
    dead_facet_owners.append( surf );
  }

PRINT_INFO("Collapsing edge (%f,%f,%f)->(%f,%f,%f) (%f)\n",
  keep_pt->coordinates().x(), keep_pt->coordinates().y(), keep_pt->coordinates().z(),
  dead_pt->coordinates().x(), dead_pt->coordinates().y(), dead_pt->coordinates().z(),
  (keep_pt->coordinates()-dead_pt->coordinates()).length());
  
    // Get PartitionCurve to update for collapsed edge.
  PartitionCurve* dead_edge_owner = 0;
  dead_edge = dead_pt->shared_edge(keep_pt);
  if (dead_edge)
    dead_edge_owner = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(dead_edge));
  
    // Collapse the edge
  if ( !keep_pt->collapse_edge(dead_pt) )
    return CUBIT_FAILURE;
  
    // Update owning curve
  if (dead_edge_owner)
    dead_edge_owner->remove_dead_facet( (CubitFacetEdgeData*)dead_edge );
  
    // Update facet owners for dead facets.
  dead_facet_ds.reset();
  dead_facet_owners.reset();
  for ( i = dead_facet_ds.size(); i--; )
  {
    CubitFacetData* facet = dead_facet_ds.get_and_step();
    PartitionSurface* facet_owner = dead_facet_owners.get_and_step();
    if (!facet_owner)
      unowned->append(facet);
    else 
      facet_owner->notify_destroyed(facet);
 }
  
  return CUBIT_SUCCESS;
}




//-------------------------------------------------------------------------
// Purpose       : Seam a list of facet edges from the boundary of a
//                 surface facetting with all the curves that are 
//                 partitions of the same real curve on the boundary of
//                 that surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/08/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::seam_curves( DLIList<PartitionCurve*>& curve_list,
                                            DLIList<CubitFacetEdge*>& edge_list,
                                            DLIList<CubitFacetData*>& facets )
{
  int i;
  if (!curve_list.size() || !edge_list.size())
    return CUBIT_FAILURE;
  
  if (DEBUG_FLAG(145))
    draw_edges( edge_list, CUBIT_WHITE, true, true, CUBIT_WHITE, false );
  
  DLIList<CubitFacetData*> old_facets, new_facets;
  DLIList<CubitFacetEdgeData*> dead_edge_ptrs;
  PartitionCurve* curve = 0;
  
    // All the curves in the attached coedge_list should belong
    // to the same real curve and thus have a single continuous 
    // parameterization over the entire set of curves.  Verify
    // this now.
  PartitionPoint *start_point, *end_point;
  Curve* real_curve = get_real_curve(curve_list, start_point, end_point);
  if (!real_curve)
  {
    PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__ );
    assert(false);
    return CUBIT_FAILURE;
  }
  
  double period;
  const bool periodic = real_curve->is_periodic(period);
  const bool fwdparam = (periodic ?  period > 0.0 : real_curve->start_param() < real_curve->end_param());
  
  edge_list.reset();
  CubitPoint *next_vtx, *curr_vtx;
  if (edge_list.size() == 1)
  {
    next_vtx = edge_list.get()->point(1);
    curr_vtx = edge_list.get()->point(0);
  }
  else
  {
    next_vtx = edge_list.get()->shared_point(edge_list.next());
    curr_vtx = edge_list.get()->other_point(next_vtx);
  }
  assert(curr_vtx&&next_vtx);
    // Special case:  If curve_list is a closed loop of curves,
    // a geometric check is required to determine if the passed
    // list of facet edges is reversed wrt the curve list.
  if (start_point == end_point)
  {
    CubitVector tangent, closest, start(curr_vtx->coordinates()), end(next_vtx->coordinates());
    real_curve->closest_point( 0.5*(start+end), closest, &tangent );
    if ((end - start) % tangent < 0.0)
      edge_list.reverse();    
  }
    // Otherwise just compare end coordinates
  else
  {
    if ((start_point->coordinates() - curr_vtx->coordinates()).length_squared() >
        (end_point->coordinates() - curr_vtx->coordinates()).length_squared())
    {
      if (edge_list.size() == 1)
      {
        next_vtx = curr_vtx;
        curr_vtx = edge_list.get()->other_point(next_vtx);
      }
      else
      {
        edge_list.reverse();
        edge_list.reset();
        next_vtx = edge_list.get()->shared_point(edge_list.next());
        curr_vtx = edge_list.get()->other_point(next_vtx);
      }
      assert(next_vtx&&curr_vtx);
    }
  }

  DLIList<CubitFacetEdgeData*> edge_data_list;
  CAST_LIST(edge_list, edge_data_list, CubitFacetEdgeData);
  assert(edge_list.size() == edge_data_list.size());
  
    // Find the subset of edges which belong on each curve.
    // Begin with the first edge in the list, and the point
    // at the end (opposite the curve start point) of the edge.
  DLIList<CubitFacetEdgeData*> curve_edge_list;  // edges for current curve
  PartitionPoint* vertex = start_point;      // end point of current curve
  curve_list.reset(); 
  edge_data_list.reset();
  edge_data_list.reverse();
  CubitFacetEdgeData* edge = edge_data_list.pop();
  CubitPoint* point = start_point->facet_point(); // edge end point
  point = edge->other_point(point);
  if (!point)
  {
    PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // Iterate over all but the last curve
  double prev_vtx_param = real_curve->start_param();
  double prev_pt_param = real_curve->start_param();
  for (i = curve_list.size() - 1; i--; )
  {
    curve_edge_list.clean_out();
    
    curve = curve_list.get_and_step();
    vertex = curve->other_point(vertex);
    const CubitVector vtx_pos(vertex->coordinates());
    double vtx_param = real_curve->u_from_position(vtx_pos);
    double point_param = real_curve->u_from_position(point->coordinates());
    if (periodic)
    {
      if ((vtx_param < prev_vtx_param) == fwdparam)
        vtx_param += period;
      if ((point_param < prev_pt_param) == fwdparam)
        point_param += period;
      prev_vtx_param = vtx_param;
      prev_pt_param = point_param;
    }
    
      // Iterate until the edge "containing" the vertex
    while ((vtx_param > point_param) == fwdparam)
    {
        // decrement edge count each time we add one of the
        // original, input edges to the cuve edge list.
      curve_edge_list.append(edge);
      
      if (!edge_data_list.size()) { assert(0); return CUBIT_FAILURE; }

        // Get the next edge and facet point
      edge = edge_data_list.pop();
      point = edge->other_point(point);
      point_param = real_curve->u_from_position(point->coordinates());
      if (periodic)
      {
        if ((point_param < prev_pt_param) == fwdparam)
          point_param += period;
        prev_pt_param = point_param;
      }
    }
  
    CubitVector edge_pos;
    edge->closest_point(vtx_pos, edge_pos);
  
    CubitFacetEdgeData* new_edge;
    CubitPoint* new_point = split_edge_closest( edge, edge_pos, 0.1*edge->length(), new_edge, facets );
    if (!new_point)
      return CUBIT_FAILURE;
    
      // If the edge was split...
    if (new_edge)
    {
        // If new_edge ends at the "end" point
        // swap the edges.
      if (new_edge->other_point(point))
        std::swap(edge, new_edge);
        // Put the earlier peice of the edge in the curve
        // list.  Keep the latter peice for the next curve
        // iteration.
      curve_edge_list.append(new_edge);
    }
      // If the input position was the same as the edge end point...
    else if (new_point == point)
    {
        // Add edge to list for current curve and decrement count.
      curve_edge_list.append(edge);
      if (!edge_data_list.size()) { assert(0); return CUBIT_FAILURE; }
      
        // Get the next edge and facet point for the next iteration.
      edge = edge_data_list.pop();
      point = edge->other_point(point);
    }
    
      // Move edge split point to location of vertex
#ifdef PART_SURF_REFACET
    double d_sqr = (vtx_pos - edge_pos).length_squared();
    bool close = d_sqr < (GEOMETRY_RESABS*GEOMETRY_RESABS);
    old_facets.clean_out();
    new_facets.clean_out();
    if (close)
      dynamic_cast<CubitPointData*>(new_point)->set(vtx_pos);
    else if (!fix_move_point( new_point, vtx_pos, facets, old_facets, new_facets ))
      return CUBIT_FAILURE;
    assert(!old_facets.size() == !new_facets.size());
    facets -= old_facets;
    facets += new_facets;
#else
    dynamic_cast<CubitPointData*>(new_point)->set(vtx_pos); 
#endif
      // Merge the points
    CubitPointData* vtx_pointd = vertex->facet_point();
    CubitPointData* new_pointd = dynamic_cast<CubitPointData*>(new_point);
    if (vtx_pointd)
      vtx_pointd->merge_points(new_pointd);
    else
      vertex->facet_point(new_pointd);
    
    if (!curve->has_facet_data())
      curve->set_facet_data(curve_edge_list);
    else if(!seam_curve( curve_edge_list, curve, facets, &dead_edge_ptrs ))
      return CUBIT_FAILURE;
      
    while (dead_edge_ptrs.size())
    {
      CubitFacetEdgeData* edge = dead_edge_ptrs.pop();
      boundary_set.erase(edge);
      interior_set.erase(edge);
    }

  } // for(curve_list)
  
    // Do last curve

  curve = curve_list.get();
  curve_edge_list.clean_out();
  if (edge)
    curve_edge_list.append(edge);
  while (edge_data_list.size())
    curve_edge_list.append(edge_data_list.pop());
  
  if(!curve_edge_list.size()) { assert(0); return CUBIT_FAILURE; }
  
  if (!curve->has_facet_data())
    curve->set_facet_data(curve_edge_list);
  else if(!seam_curve( curve_edge_list, curve, facets, &dead_edge_ptrs ))
    return CUBIT_FAILURE;
      
  while (dead_edge_ptrs.size())
  {
    CubitFacetEdgeData* edge = dead_edge_ptrs.pop();
    boundary_set.erase(edge);
    interior_set.erase(edge);
  }
    
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Seam a new list of facet edges with the list of facet
//                 edges on a curve.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::seam_curve( DLIList<CubitFacetEdgeData*>& edge_list,
                                           PartitionCurve* curve,
                                           DLIList<CubitFacetData*>& facets,
                                           DLIList<CubitFacetEdgeData*>* dead_ptrs )
{
  DLIList<CubitFacetData*> old_facets, new_facets;
  const double FRACT_TOL = 0.1;

  int curve_color = CUBIT_GREEN;
  if (DEBUG_FLAG(145))
  {
    RefEntity* ent = dynamic_cast<RefEntity*>(curve->topology_entity());
    while (ent && ent->color() < 0)
    {
      DLIList<RefEntity*> list;
      ent->get_parent_ref_entities(list);
      list.reset();
      if (!list.size()) break;
      ent = list.get();
    }
    if (ent)
      curve_color = ent->color();
  }
  
    // Vertices should already be "seamed".  Verify...
    
  PartitionPoint* start_vtx = curve->start_point();
  PartitionPoint* end_vtx = curve->end_point();
  CubitPointData* start_point = start_vtx->facet_point();
  CubitPointData* end_point = end_vtx->facet_point();
  
  edge_list.last();
  CubitFacetEdgeData* last_edge = edge_list.get();
  edge_list.reset();
  CubitFacetEdgeData* first_edge = edge_list.get();
  if (!first_edge->other_point(start_point) || 
      !last_edge->other_point(end_point))
  {
    PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // Get direction of parameter (probably always increasing...)
  double period;
  const bool fwdparam = curve->is_periodic(period) ? period > 0.0 :
                        curve->start_param() < curve->end_param();
  
    // Get list of facet edges on curve.
  DLIList<CubitFacetEdgeData*> curve_edges;
  curve->get_facet_data(curve_edges);
  if (!curve_edges.size())
  {
    curve->set_facet_data( edge_list );
    return CUBIT_SUCCESS;
  }
  
    // Seam edge lists
  double curve_param = curve->start_param();
  double new_param = curve->start_param();
  curve_edges.reset();
  edge_list.reset();
  CubitFacetEdgeData* curve_edge = curve_edges.get_and_step();
  CubitFacetEdgeData* new_edge = edge_list.get_and_step();
  CubitPoint* prev_point = start_point;
  double orig_crv_len = curve_edge->length();
  double orig_new_len = new_edge->length();
  double split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
  while (curve_edge || new_edge)
  {
    if (!curve_edge || !new_edge)
    {
      PRINT_ERROR("Internal error at %s:%d\n", __FILE__, __LINE__ );
      assert(0);
      return CUBIT_FAILURE;
    }
  
    CubitPoint* curve_point = curve_edge->other_point(prev_point);
    CubitPoint* new_point = new_edge->other_point(prev_point);
    if (DEBUG_FLAG(145)) {
      GfxDebug::draw_point(curve_point->coordinates(), curve_color);
      GfxDebug::draw_facet_edge(curve_edge, curve_color);
      GfxDebug::draw_point(new_point->coordinates(), curve_color + 1);
      GfxDebug::draw_facet_edge(new_edge, curve_color + 1);
      GfxDebug::flush();
    }
    
    double prev_curve_param = curve_param;
    double prev_new_param = new_param;
    curve_param = curve->u_from_position(curve_point->coordinates());
    new_param = curve->u_from_position(new_point->coordinates());
    if (curve->is_periodic(period))
    {
      if ((prev_curve_param > curve_param) == fwdparam)
        curve_param += period;
      if ((prev_new_param > new_param) == fwdparam)
        new_param += period;
    }
    
    if ((curve_param > new_param) == fwdparam)
    {
      CubitFacetEdgeData* split_edge = 0;
      CubitVector edge_pos;
      curve_edge->closest_point( new_point->coordinates(), edge_pos );
      CubitPoint* split = split_edge_closest( curve_edge, 
                                              edge_pos, 
                                              split_tol,
                                              split_edge, 
                                              facets );
      if(!split)
        return CUBIT_FAILURE;
      
      if (split == prev_point)
      {
          // new_edge is very small compared to split_edge
        assert(!split_edge);
        old_facets.clean_out();
        if (!collapse_edge(prev_point, new_point, &old_facets))
        {
          if (!collapse_edge(new_point, prev_point, &old_facets))
          {
            assert(0);
            return CUBIT_FAILURE;
          }
          std::swap(prev_point, new_point);
        }
        if (dead_ptrs)
          dead_ptrs->append( new_edge );
        facets -= old_facets;
      
        if (edge_list.is_at_beginning())
          new_edge = 0;
        else
        {
          new_edge = edge_list.get_and_step();
          orig_new_len = new_edge->length();
          split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
        }
         
        continue;
      }
      
      if (split_edge)
      {
        if (curve_edge->other_point(prev_point))
          std::swap(split_edge, curve_edge);
      }
      else
      {
        assert(split == curve_point);
        split_edge = curve_edge;
        if (curve_edges.is_at_beginning())
          curve_edge = 0;
        else
        {
          curve_edge = curve_edges.get_and_step();
          orig_crv_len = curve_edge->length();
          split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
        }
      }
      
        // Move split point to position of other point (safely)
#ifdef PART_SURF_REFACET
      double d_sqr = (edge_pos - new_point->coordinates()).length_squared();
      bool close = d_sqr < (GEOMETRY_RESABS*GEOMETRY_RESABS);
      old_facets.clean_out();
      new_facets.clean_out();
      if (close)
        dynamic_cast<CubitPointData*>(split)->set(new_point->coordinates());
      else if (!fix_move_point( split, new_point->coordinates(), facets, old_facets, new_facets ))
        return CUBIT_FAILURE;
      facets -= old_facets;
      facets += new_facets;
#else
      dynamic_cast<CubitPointData*>(split)->set(new_point->coordinates());
#endif
      
      CubitStatus s1 = split->merge_points(new_point);
      CubitStatus s2 = split_edge->merge_edges(new_edge);
      if (dead_ptrs)
        dead_ptrs->append(new_edge);
      assert(s1 && s2);
      prev_point = split;
      
      if (edge_list.is_at_beginning())
        new_edge = 0;
      else
      {
        new_edge = edge_list.get_and_step();
        orig_new_len = new_edge->length();
        split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
      }  
    }
    else
    {
      CubitFacetEdgeData* split_edge = 0;
      CubitVector edge_pos;
      new_edge->closest_point( curve_point->coordinates(), edge_pos );
      CubitPoint* split = split_edge_closest( new_edge, 
                                              edge_pos, 
                                              split_tol,
                                              split_edge,
                                              facets );
      if(!split)
        return CUBIT_FAILURE;
      
      if (split == prev_point)
      {
          // curve_edge is very small compared to split_edge
        assert(!split_edge);
        old_facets.clean_out();
        if (!collapse_edge(prev_point,curve_point,&old_facets))
        {
          if (!collapse_edge(curve_point, prev_point, &old_facets))
          {
            assert(0);
            return CUBIT_FAILURE;
          }
          std::swap(curve_point, prev_point);
        }
        facets -= old_facets;
        if (dead_ptrs)
          dead_ptrs->append(curve_edge);
      
        if (curve_edges.is_at_beginning())
          curve_edge = 0;
        else
        {
          curve_edge = curve_edges.get_and_step();
          orig_crv_len = curve_edge->length();
          split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
        }  
        continue;
      }
        
      if (split_edge)
      {
        if (new_edge->other_point(prev_point))
          std::swap(split_edge, new_edge);
      }
      else
      {
        assert(split == new_point);
        split_edge = new_edge;
        if (edge_list.is_at_beginning())
          new_edge = 0;
        else
        {
          new_edge = edge_list.get_and_step();
          orig_new_len = new_edge->length();
          split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
        }  
      }    
      
        // Move split point to position of other point (safely)
#ifdef PART_SURF_REFACET
      double d_sqr = (edge_pos - curve_point->coordinates()).length_squared();
      bool close = d_sqr < (GEOMETRY_RESABS*GEOMETRY_RESABS);
      old_facets.clean_out();
      new_facets.clean_out();
      if (close)
        dynamic_cast<CubitPointData*>(split)->set(curve_point->coordinates());
      else if (!fix_move_point( split, curve_point->coordinates(), facets, old_facets, new_facets ))
        return CUBIT_FAILURE;
      facets -= old_facets;
      facets += new_facets;
#else
      dynamic_cast<CubitPointData*>(split)->set(curve_point->coordinates());
#endif
      
      CubitStatus s2, s1 = CUBIT_SUCCESS;
      if (curve_point != split)
        s1 = curve_point->merge_points(split);
      s2 = curve_edge->merge_edges(split_edge);
      if (dead_ptrs)
        dead_ptrs->append(split_edge);
      assert(s1 && s2);
      prev_point = curve_point;
      
      if (curve_edges.is_at_beginning())
        curve_edge = 0;
      else
      {
        curve_edge = curve_edges.get_and_step();
        orig_crv_len = curve_edge->length();
        split_tol = FRACT_TOL * CUBIT_MIN(orig_crv_len, orig_new_len);
      }  
    }
  }
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Retriangulate if necessary for moving a boundary point
//                 into the interior of a facet patch.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/07/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::fix_move_point( 
                                      CubitPoint* point,
                                      const CubitVector& new_pos,
                                      const DLIList<CubitFacetData*>& facetds,
                                      DLIList<CubitFacetData*>& old_facets,
                                      DLIList<CubitFacetData*>& new_facets,
                                      PartitionSurface* surface )
{
  int i;
  DLIList<CubitFacet*> facets(facetds.size());
  for (i = 0; i < facetds.size(); i++)
    facets.append(facetds.next(i));
  
    // Look for a triangle containing the passed point.
  CubitVector closest_pos, vect;
  CubitFacet* facet = closest_facet(new_pos, facetds, closest_pos );
  if ( facet == NULL)
    return CUBIT_FAILURE;

  if (facet->point_index(point) >= 0)
  {
    dynamic_cast<CubitPointData*>(point)->set(new_pos);
    return CUBIT_SUCCESS;
  }

  PRINT_WARNING("Refacetting surface near (%f,%f,%f)\n",
    new_pos.x(), new_pos.y(), new_pos.z());

    // Line direction
  const CubitVector dir(point->coordinates() - new_pos);

  if (DEBUG_FLAG(145)) 
  { 
    GfxDebug::draw_facet(facet, CUBIT_BLUE); 
    GfxDebug::draw_point(point->coordinates(), CUBIT_BLUE);
    GfxDebug::draw_point(new_pos, CUBIT_CYAN);
    GfxDebug::draw_line(point->coordinates(), new_pos, CUBIT_BLUE);
    GfxDebug::flush(); 
  }
  
    // Find where line exits first triangle
  CubitFacetEdge* edge = 0;
  bool leftofpoint[3];
  for ( i = 0; i < 3; i++)
  {
    CubitVector vect = facet->point(i)->coordinates() - new_pos;
    vect *= dir;
    leftofpoint[i] = (vect % facet->normal()) >= 0.0;
  }
  
  for ( i = 0; i < 3; ++i)
    if (leftofpoint[i] && !leftofpoint[(i+1)%3])
    {
      edge = facet->edge((i+2)%3);
      break;
    }
  
  if (!edge)
    return CUBIT_FAILURE;
  
    // Find the set of triangles the line from the old 
    // position to the new position crosses.
  DLIList<CubitFacet*> intersect_list, edge_facets;
  while (facet->point_index(point) < 0)
  {
    intersect_list.append(facet);

    if (DEBUG_FLAG(145)) 
    { 
      GfxDebug::draw_facet_edge(edge, CUBIT_WHITE); 
      GfxDebug::draw_facet(facet, CUBIT_BLUE); 
      GfxDebug::flush(); 
    }
      
      // Find the facet on the other side of the edge
    edge_facets.clean_out();
    if (surface)
      PartSurfFacetTool::edge_facets(surface, edge, edge_facets);
    else
      PartSurfFacetTool::edge_facets( edge, facets, edge_facets );
    edge_facets.move_to(facet);
    assert(edge_facets.get() == facet);

    if(edge_facets.size() != 2)
      break;
    
    facet = edge_facets.step_and_get();
    
      // Find the edge that the line intersects exiting the facet
    i = facet->edge_index(edge);
    vect = facet->point(i)->coordinates() - new_pos;
    vect *= dir;
    if (vect % facet->normal() >= 0.0)
      edge = facet->edge((i+1)%3);
    else
      edge = facet->edge((i+2)%3);
  }
  
  if (facet)
  {
    if (DEBUG_FLAG(145)) 
    { 
      GfxDebug::draw_facet_edge(edge, CUBIT_WHITE); 
      GfxDebug::flush(); 
    }
    intersect_list.append_unique(facet);
  }
  
    // Get ordered boundary
  DLIList<CubitPoint*> boundary_list;
  DLIList<CubitFacetEdge*> point_edges;
  CubitFacetEdge* prev_edge = 0;
  CubitFacet* prev_facet = 0;
    // Find one boundary edge to start with
  for (i = intersect_list.size(); !prev_edge && i--; )
  {
    CubitFacet* facet = intersect_list.get_and_step();
    for (int j = 0; j < 3; j++)
    {
      CubitFacetEdge* edge = facet->edge(j);
      edge_facets.clean_out();
     if (surface)
       PartSurfFacetTool::edge_facets( surface, edge, edge_facets);
     else
       PartSurfFacetTool::edge_facets( edge, facets, edge_facets );
      edge_facets.intersect(intersect_list);
      if (edge_facets.size() == 1)
      {
        prev_edge = edge;
        prev_facet = facet;
        break;
      }
    }
  }
  CubitFacetEdge* first_edge = prev_edge;
    // Traverse boundary edges in order
  do 
  {
    int direction = prev_facet->edge_use( prev_facet->edge_index(prev_edge) );
    CubitPoint* pt = prev_edge->point( direction == -1 ? 0 : 1 );
    point_edges.clean_out();
    pt->edges( point_edges );
    
    if (DEBUG_FLAG(145))
    {
      GfxDebug::draw_facet_edge( prev_edge, CUBIT_LIGHTBLUE );
      GfxDebug::flush();
    }
    boundary_list.append(pt);
    
    point_edges.move_to(prev_edge);
    assert(point_edges.get() == prev_edge);
    point_edges.extract();
    prev_edge = 0;

    while (point_edges.size())
    {
      CubitFacetEdge* edge = point_edges.pop();
      edge_facets.clean_out();
      if (surface)
        PartSurfFacetTool::edge_facets( surface, edge, edge_facets);
      else
        PartSurfFacetTool::edge_facets( edge, facets, edge_facets );
      edge_facets.intersect(intersect_list);
      if (edge_facets.size() == 1)
      {
        if (prev_edge) {
          prev_edge = 0;
          break;
        }
        prev_edge = edge;
        prev_facet = edge_facets.get();
      }
    }
  
    if (!prev_edge)
    {
      PRINT_ERROR("Problems finding boundary to re-facet around split point.\n");
      return CUBIT_FAILURE;
    }
  } while (prev_edge != first_edge);

 
    // Construct new facets in affected region by connecting 
    // the split point to each boundary point
  int junk = 0;
  //DLIList<CubitFacetData*> new_facets, old_facets(intersect_list.size());
  if (!boundary_list.move_to(point))
  {
    PRINT_ERROR("Problems finding region around split point for re-facetting.\n");
    return CUBIT_FAILURE;
  }
  boundary_list.step();
  CubitPoint* prev_pt = boundary_list.get_and_step();
  for (i = boundary_list.size(); i > 2; i--)
  {
    CubitPoint* next_pt = boundary_list.get_and_step();
    CubitFacetData* new_facet = new CubitFacetData( point, prev_pt, next_pt, &junk);
    new_facets.append(new_facet);
    prev_pt = next_pt;
    
    if (DEBUG_FLAG(145))
    {
      GfxDebug::draw_facet( new_facet, CUBIT_YELLOW );
      GfxDebug::flush();
    }
  }
  
  CAST_LIST( intersect_list, old_facets, CubitFacetData );
  assert(intersect_list.size() == old_facets.size());
  dynamic_cast<CubitPointData*>(point)->set(new_pos);
  //replace_facets( old_facets, new_facets );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Find closest facet and point on facet
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/28/03
//-------------------------------------------------------------------------
CubitFacet* PartSurfFacetTool::closest_facet(
                                   const CubitVector& input_position,
                                   const DLIList<CubitFacetData*>& facet_list,
                                   CubitVector& result_position )
{
  CubitFacet* return_val = 0;
  double closest_dist_sqr = CUBIT_DBL_MAX;
  CubitVector facet_position;
  
  const int size = facet_list.size();
  for ( int i = 0; i < size; i++ ) {
    CubitFacet* facet = facet_list.next(i);
    closest_pt_on_facet( facet, input_position, facet_position );
    double dist_sqr = (input_position - facet_position).length_squared();
    if ( dist_sqr < closest_dist_sqr ) {
      closest_dist_sqr = dist_sqr;
      result_position = facet_position;
      return_val = facet;
    }
  }
  
  return return_val;
}

//-------------------------------------------------------------------------
// Purpose       : Return edge start/end if input position is sufficiently
//                 close.  Otherwise split the edge and return the new
//                 point.
//
// Special Notes : For all facets modified by this operation, owning surfaces
//                 will be notified of the change.  For any un-owned facets,
//                 new facets will be appended to the passed list.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/01/03
//-------------------------------------------------------------------------
CubitPoint* PartSurfFacetTool::split_edge_closest( 
                                     CubitFacetEdgeData* old_edge,
                                     const CubitVector& pos,
                                     double tolerance,
                                     CubitFacetEdgeData*& new_edge,
                                     DLIList<CubitFacetData*>& new_facets )
{
    // square of absolute tolerance
  const double GEOM_TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS; 
    // tolerance as fraction of edge length.
  const double TOL_SQR = tolerance < GEOM_TOL_SQR ? GEOM_TOL_SQR : tolerance * tolerance;
  
  new_edge = 0;

  const CubitVector start(old_edge->point(0)->coordinates());
  const CubitVector end(old_edge->point(1)->coordinates());
//   const CubitVector dir(end - start);
  //const double fract = (dir % (pos - start)) / dir.length_squared();

    // If near the start/end of the edge and start/end point is
    // not already owned by some other vertex
  const double start_sqr = (start - pos).length_squared();
  const double end_sqr   = (end   - pos).length_squared();
  if (start_sqr < TOL_SQR && !TDVGFacetOwner::get(old_edge->point(0)))
    return old_edge->point(0);
  else if (end_sqr < TOL_SQR && !TDVGFacetOwner::get(old_edge->point(1)))
    return old_edge->point(1);

    // Else if start/end of edge are same position as vertex
  else if(start_sqr < GEOM_TOL_SQR)
    return old_edge->point(0);
  else if(end_sqr < GEOM_TOL_SQR)
    return old_edge->point(1);

    // Otherwise split the edge
  CubitFacet* a_face = old_edge->adj_facet(0);
  CubitFacetData* facet = dynamic_cast<CubitFacetData*>(a_face);
  if (!facet)
    return 0;

  //CubitVector edge_pos((1.0-edge_param)*edge_start + edge_param*edge_end);
  CubitPoint* start_pt = old_edge->point(0);
  CubitPoint*   end_pt = old_edge->point(1);
  CubitPoint*   new_pt = facet->split_edge( start_pt, end_pt, pos );
  CubitFacetEdge* temp_edge = new_pt->shared_edge(end_pt);
  new_edge = dynamic_cast<CubitFacetEdgeData*>(temp_edge);
  if (!new_edge)
    return 0;
  
    // Update partition curve facetting
  PartitionCurve* edge_owner = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(old_edge));
  if (edge_owner)
    edge_owner->notify_split(old_edge, new_edge);

    // Update partition surfaces for split of adjacent facets
  int i, junk = -1;
  DLIList<CubitFacet*> facets(old_edge->num_adj_facets());
  old_edge->facets(facets);
  
  for ( i = facets.size(); i--;  )
  {
    a_face = facets.get_and_step();
    facet = dynamic_cast<CubitFacetData*>(a_face);
    assert(!!facet);
    int edge_index = facet->edge_index( start_pt, new_pt, junk );
    CubitPoint* other_pt = facet->point( edge_index );
    if (!other_pt->shared_edge( new_pt ))
      new CubitFacetEdgeData( other_pt, new_pt );

    CubitFacet* a_face = facet->shared_facet( new_pt, other_pt );
    CubitFacetData* other_facet = dynamic_cast<CubitFacetData*>(a_face);
    if (!other_facet) { assert(0); continue; }

    PartitionEntity* owner = TDVGFacetOwner::get(facet);
    PartitionSurface* surf = dynamic_cast<PartitionSurface*>(owner);
    assert(!owner || surf);
    
    if (owner)
      owner->notify_split( facet, other_facet );
    else
      new_facets.append(other_facet);
  }
  
  return new_pt;
}

//-------------------------------------------------------------------------
// Purpose       : Spatially match a set of geometric points one-to-one
//                 with a subset of a set of facet vertices.
//
// Special Notes : If the geometric point already has an attached facet
//                 point, the points will be merged such that the facet
//                 point originally on the geometric point survives the
//                 merge.  If the geometric point does not already have
//                 a facet point, the input facet point will be attached
//                 to the geometric point.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/08/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::associate_points( 
                                     DLIList<CubitPoint*>& facet_points,
                                     DLIList<PartitionPoint*>& geom_points )
{
  int i;
  DLIList<int> geom_indices(geom_points.size());
  
  for (i = 0; i < geom_points.size(); i++ )
    geom_indices.append(i);

  for (i = facet_points.size(); i--; )
    facet_points.step_and_get()->marked(-1);

  geom_points.reset();
  facet_points.reset();
  while (geom_indices.size())
  {
    int index = geom_indices.pop();
    PartitionPoint* geom_pt = geom_points.next(index);
    const CubitVector pos(geom_pt->coordinates());
    double shortest_dist_sqr = CUBIT_DBL_MAX;
    int closest_index = -1;
    for (i = 0; i < facet_points.size(); i++ )
    {
      CubitPoint* facet_pt = facet_points.next(i);
      double dist_sqr = (pos - facet_pt->coordinates()).length_squared();
      if (dist_sqr < shortest_dist_sqr)
      {
        if (facet_pt->marked() == -1)
        {
          shortest_dist_sqr = dist_sqr;
          closest_index = i;
        }
        else 
        {
          int j = facet_pt->marked();
          CubitPoint* other_point = facet_points.next(j);
          double other_dist_sqr = (pos - other_point->coordinates()).length_squared();
          if (other_dist_sqr > dist_sqr)
          {
            shortest_dist_sqr = dist_sqr;
            closest_index = i;
            facet_pt->marked(-1);
            geom_indices.append(j);
          }
        }
      }
    }
    
    if (closest_index < 0)
    {
      assert(false);
      for (i = facet_points.size(); i--; )
        facet_points.step_and_get()->marked(0);
      for (i = geom_points.size(); i--; )
        geom_points.step_and_get()->mark = 0;
      return CUBIT_FAILURE;
    }
    geom_pt->mark = closest_index;
    facet_points.next(closest_index)->marked(index);
  }
  
  for (i = facet_points.size(); i--; )
    facet_points.step_and_get()->marked(0);
  
  facet_points.reset();
  for (i = geom_points.size(); i--; )
  {
    PartitionPoint* geom_pt = geom_points.step_and_get();
    CubitPoint* facet_pt = facet_points.next(geom_pt->mark);
    geom_pt->mark = 0;
    
    assert(!TDVGFacetOwner::get(facet_pt));
    facet_pt->set(geom_pt->coordinates());
    
    CubitPointData* point_data = dynamic_cast<CubitPointData*>(facet_pt);
    assert(!!point_data);
    
    if (geom_pt->facet_point())
      geom_pt->facet_point()->merge_points( point_data );
    else
      geom_pt->facet_point(point_data);
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get points and edges in facet patch, separated into
//                 boundary and internal sets where a boundary point is
//                 a point for which one or more adjacent edges is on the
//                 boundary of the patch.
//
// Special Notes : Aborts and returns failure if the input facets are not
//                 a manifold patch.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/08/03
//-------------------------------------------------------------------------
CubitStatus PartSurfFacetTool::get_facet_points_and_edges(
                                   const DLIList<CubitFacetData*>& facets,
                                   DLIList<CubitPoint*>& boundary_points,
                                   DLIList<CubitPoint*>& interior_points,
                                   DLIList<CubitFacetEdge*>& boundary_edges,
                                   DLIList<CubitFacetEdge*>& interior_edges )
{
  int i, j;
  CubitStatus result = CUBIT_SUCCESS;
  CubitFacetData* facet;
  DLIList<CubitFacetEdge*> edge_list;
  const int count = facets.size();
  
    // Initialize marks
  for (i = 0; i < count; i++ )
  {
    facet = facets.next(i);
    for (j = 0; j < 3; j++)
    {
      facet->edge(j)->marked(0);
      facet->point(j)->marked(1);
    }
  }
  
    // Mark each edge with a count of the number of
    // adjcent facets.
  for (i = 0; i < count; i++ )
  {
    facet = facets.next(i);
    for (j = 0; j < 3; j++)
      facet->edge(j)->marked(facet->edge(j)->marked() + 1);
  }
  
    // Collect points and clear point marks
  for (i = 0; i < count; i++ )
  {
    facet = facets.next(i);
    for (j = 0; j < 3; j++)
    {
      CubitPoint* point = facet->point(j);
      if (!point->marked())
        continue;
      point->marked(0);
      
      point->edges(edge_list);
      bool boundary = false;
      while (edge_list.size())
        if (edge_list.pop()->marked() == 1)
          boundary = true;
      if (boundary)
        boundary_points.append(point);
      else 
        interior_points.append(point);
    }
  }
  
    // Collect edges and clear edge marks   
  for (i = 0; i < count; i++ )
  {
    facet = facets.next(i);
    for (j = 0; j < 3; j++)
    {
      CubitFacetEdge* edge = facet->edge(j);
      switch (edge->marked())
      {
        case 0:
          continue;
        case 1:
          boundary_edges.append(edge);
          break;
        case 2:
          interior_edges.append(edge);
          break;
        default:
          result = CUBIT_FAILURE;
          assert(false);
      }
      edge->marked(0);
    }
  }
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Get facets adjacent to edge and owned by this surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/07/03
//-------------------------------------------------------------------------
void PartSurfFacetTool::edge_facets( PartitionSurface* surf,
                                     CubitFacetEdge* edge,
                                     DLIList<CubitFacet*>& facets )
{
  edge->facets(facets);

  for (int i = facets.size(); i--; )
    if (TDVGFacetOwner::get(facets.step_and_get()) != surf)
      facets.change_to(0);
  
  facets.remove_all_with_value(0);
}



//-------------------------------------------------------------------------
// Purpose       : Get facets in passed list adjacent to passed edge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/07/03
//-------------------------------------------------------------------------
void PartSurfFacetTool::edge_facets( CubitFacetEdge* edge,
                                    const DLIList<CubitFacet*>& input_facets,
                                    DLIList<CubitFacet*>& output_facets )
{
  edge->facets(output_facets);
  output_facets.intersect( input_facets );
}
