//-------------------------------------------------------------------------
// Filename      : CompSurfFacets.cpp
//
// Purpose       : Encapsulate facet date used to speed up composite surfaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/20/03
//-------------------------------------------------------------------------

#include "CompSurfFacets.hpp"
#include "Surface.hpp"
#include "GeometryQueryEngine.hpp"
#include "GMem.hpp"
#include "OctTree.hpp"
#include "GfxDebug.hpp"

CubitStatus CompSurfFacets::setup( const SurfPtrList& surface_data )
{
  GMem gmem;
  int num_points, num_facets, num_triangles;
  
  // Set the size of this array of flags.
  ignoreFlags.resize(surface_data.size());

  for ( unsigned int index = 0; index < surface_data.size(); index++ )
  {
    ignoreFlags[index] = 0;
    Surface* surface = surface_data[index];
    
    CubitStatus result = surface->get_geometry_query_engine()
      ->get_graphics( surface, num_triangles, num_points, num_facets, &gmem );
    if ( !result )
      return CUBIT_FAILURE;
      
    int* f_itor = gmem.facet_list();
    int* f_end = f_itor + num_facets;
    int tri_count = 0;
    while ( f_itor < f_end )
    {
        // If facet is not a triangle, it will be split into
        // triangles later.  Update tri_count accordingly.
      int pt_count = *f_itor;
      tri_count += pt_count - 2;
      f_itor += (pt_count + 1);
    }
    
    if (tri_count != num_triangles)
      PRINT_WARNING("Non-triangular facets encountered in surface graphics.\n");
    
    assert(gmem.fListCount == num_facets);
    assert(gmem.pointListCount == num_points);
    
    int offset = pointList.size();
    pointList.resize( offset + num_points );
    
    int tri_offset = triangleData.size();
    triangleData.resize( tri_offset + tri_count * 4 );
    
    PointList::iterator p_itor = pointList.begin() + offset;
    GPoint* g_itor = gmem.point_list();
    GPoint* g_end = g_itor + num_points;
    for ( ; g_itor != g_end; ++g_itor )
    {
      //CubitVector vect(g_itor->x, g_itor->y, g_itor->z);
      //CubitVector closest;
      //surface->closest_point( vect, &closest );
      //*p_itor = closest;
      p_itor->set(g_itor->x, g_itor->y, g_itor->z);
      ++p_itor;
    }
    
    IntegerList::iterator t_itor = triangleData.begin() + tri_offset;
    f_itor = gmem.facet_list();
    f_end = f_itor + num_facets;
    while ( f_itor != f_end )
    {
        // Copy facet into triangleData.
        // If facet is non-triangular, assume facet is convex
        // and reduce it to triangles.
      int count = *f_itor++;
      int base = *f_itor++;
      int prev = *f_itor++;
      for (int i = 0; i < count - 2; i++ )
      {
        *t_itor++ = index;
        *t_itor++ = offset + base;
        *t_itor++ = offset + prev;
        prev = *f_itor++;
        *t_itor++ = offset + prev;
      }
    }
  }
  
  pointsConsolidated = false;
  return CUBIT_SUCCESS;
}

void CompSurfFacets::set_ignore_flag(int index, int flag)
{
   ignoreFlags[index] = flag;
}

int CompSurfFacets::get_ignore_flag(int index)
{
   return ignoreFlags[index];
}

int CompSurfFacets::closest_index( const CubitVector& pos,
                                   CubitVector* closest_pos ) const
{
  IntegerList::const_iterator itor = triangleData.begin();
  IntegerList::const_iterator end  = triangleData.end();
  double closest_dist = CUBIT_DBL_MAX;
  double dist = CUBIT_DBL_MAX;
  int closest_index = -1;
  CubitVector close_pt;
  bool interior = false;
  
  while ( itor != end )
  {
    // If the parent CompositeSurface has been told to
    // ignore any of its surfaces during evaluation we 
    // handle it here.  If we are evaluating facets
    // that are for a surface to ignore just set the
    // distance to be the current closest distance 
    // and continue on.
    if(!ignoreFlags[*itor])
    {
      closest_tri_point( itor, pos, close_pt, interior );
      dist = (pos - close_pt).length_squared();
    }
    else
    {
      dist = closest_dist;
    }
    if ( dist < closest_dist )
    {
      closest_dist = dist;
      closest_index = *itor;
      if ( closest_pos )
        *closest_pos = close_pt;
        
      if ( dist < (GEOMETRY_RESABS*GEOMETRY_RESABS) )
        break;
    }
    itor += 4;
  }

  return closest_index;
}

int CompSurfFacets::closest_index( const CubitVector& pos,
                                   DLIList<int>& index_list,
                                   CubitVector* closest_pos ) const
{
  IntegerList::const_iterator itor = triangleData.begin();
  IntegerList::const_iterator end  = triangleData.end();
  double closest_dist = CUBIT_DBL_MAX;
  double dist = CUBIT_DBL_MAX;
  int closest_index = -1;
  CubitVector close_pt;
  bool interior = false;
  
  index_list.clean_out();
  
  while ( itor != end )
  {
    // If the parent CompositeSurface has been told to
    // ignore any of its surfaces during evaluation we 
    // handle it here.  If we are evaluating facets
    // that are for a surface to ignore just set the
    // distance to be the current closest distance 
    // and continue on.
    if(!ignoreFlags[*itor])
    {
      closest_tri_point( itor, pos, close_pt, interior );
      dist = (pos - close_pt).length_squared();
    }
    else
    {
      dist = closest_dist;
    }
    if ( dist < closest_dist )
    {
      closest_dist = dist;
      closest_index = *itor;
      if ( closest_pos )
        *closest_pos = close_pt;

      if (interior && dist < (GEOMETRY_RESABS*GEOMETRY_RESABS) )
      {
        index_list.clean_out();
        break;
      }
    }
    if ( dist < (GEOMETRY_RESABS*GEOMETRY_RESABS) ) 
      index_list.append( *itor );

    itor += 4;
  }

    // didn't find exact match
  if ( index_list.size() == 0 )
    index_list.append( closest_index );
  else
    index_list.uniquify_unordered();
  
  return index_list.size();
}

void CompSurfFacets::closest_tri_point( IntegerList::const_iterator facet,
                                        const CubitVector& p,
                                        CubitVector& result,
                                        bool& interior ) const
{
    // Get triangle vertices
    // First value in vector is the surface index.
    // Skip the surface index and get the positions
    // for the three corner indices.
  CubitVector p0(pointList[*++facet]);
  CubitVector p1(pointList[*++facet]);
  CubitVector p2(pointList[*++facet]);

//   if ( DEBUG_FLAG(87) )
//   {
//     GfxDebug::display_all();
//     GfxDebug::draw_point(p,CUBIT_RED);
//     GfxDebug::draw_line( p0, p1, CUBIT_BLUE );
//     GfxDebug::draw_line( p0, p2, CUBIT_BLUE );
//     GfxDebug::draw_line( p2, p1, CUBIT_BLUE );
//     GfxDebug::flush();
//   }

//=============================================================================
//   Algorithm from:
//     "Distance Between Point and Triangle in 3D"
//     David Eberly
//     Magic Software, Inc.
//     Sept. 28, 1999
//
//   Use barycentric coordinates.  Coordinates are
//   calculated in the range [0,det] rather than [0,1]
//   to avoid the fp division entirely where it can
//   be avoided.
//
//     ^v*t                                 
//  \R2|                 
//   \ |                 
//    \|                 
//     *p2               
//     |\                
//     | \     R1        
//  R3 |  \              
//     |   \             
//     | R0 \            
//     |     \p1         
//  ---*------*--->u*s   
//     |p0     \  R6     
//  R4 |   R5   \        
//     |         \u+v=det  
//
//=============================================================================

  CubitVector s(p0, p1); // the u (or s) axis in parameterized space
  CubitVector t(p0, p2); // the v (or t) axis in parameterized space
  CubitVector d(p,  p0); // vector from input position to corner at (u,v)
                         // = (0,0)

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

    // is interior to the boundary
  interior = false;
  if (u > GEOMETRY_RESABS && v > GEOMETRY_RESABS &&
      u+v < det - GEOMETRY_RESABS)
    interior = true;

//   if ( DEBUG_FLAG(87) )
//   {
//     GfxDebug::draw_vector( p, result, CUBIT_RED );
//     GfxDebug::flush();
//   }
}


void CompSurfFacets::debug_draw_facets() const
{
  IntegerList::const_iterator itor = triangleData.begin();
  IntegerList::const_iterator end  = triangleData.end();
  
  while ( itor != end )
  {
    int color = (1 + *itor++) % 100;
    CubitVector p0(pointList[*itor++]);
    CubitVector p1(pointList[*itor++]);
    CubitVector p2(pointList[*itor++]);
    GfxDebug::draw_line( p0, p1, color );
    GfxDebug::draw_line( p1, p2, color );
    GfxDebug::draw_line( p2, p0, color );
  }
  GfxDebug::flush();
}

void CompSurfFacets::consolidate_points( double tolerance )
{
  if (pointList.size() < 1000)
    consolidate_few_points(tolerance);
  else
    consolidate_many_points(tolerance);
}

void CompSurfFacets::consolidate_few_points( double tolerance )
{
  const double tolsqr = tolerance * tolerance;
  int* index_map = new int[pointList.size()];
    
    // Consolidate the point list.  index_map is used
    // to maintain a map between the old index of a point
    // (the index into index_map) and the new index of a
    // point (the value in index_map).
  int write = 0, read, comp, size = pointList.size();
  for( read = 0; read < size; read++ )
  {
    const CubitVector& pti = pointList[read];
    for( comp = 0; comp < write; comp++ )
    {
      const CubitVector& ptj = pointList[comp];
      if( (pti-ptj).length_squared() <= tolsqr )
        break;
    }
    
    index_map[read] = comp;
    if( comp == write )
    {
      pointList[comp] = pointList[read];
      write++;
    }
  }
  pointList.resize(write);  

    // Update the facet list using values from index_map.
  IntegerList::iterator itor = triangleData.begin();
  IntegerList::iterator end = triangleData.end();
  while( itor++ != end )
    for( int count = 3; count--; itor++ )
      *itor = index_map[*itor];

  delete [] index_map;
  pointsConsolidated = true;
}

class CubitVectOctTreeEval {
public:
static inline const CubitVector& coordinates(CubitVector* p)
  { return *p; };
};

void CompSurfFacets::consolidate_many_points( double tolerance )
{
  const double tolsqr = tolerance*tolerance;
  
    // Copy positions into array and build OctTree
  DLIList<CubitVector*> point_list(pointList.size());
  CubitVector* point_array = new CubitVector[pointList.size()];
  PointList::iterator p_itor = pointList.begin();
  PointList::iterator p_end  = pointList.end();
  CubitVector* a_itor = point_array;
  while (p_itor != p_end)
  {
    *a_itor = *p_itor;
    point_list.append(a_itor);
    ++a_itor;
    ++p_itor;
  }
  OctTree<CubitVector,CubitVectOctTreeEval> tree( point_list, tolerance );
  
    // Consolidate the point list.  index_map is used
    // to maintain map between old and new point indices.
  int* index_map = new int[pointList.size()];
  int read, write = 0, size = pointList.size();
  for (read = 0; read < size; read++)
  {
    index_map[read] = write;
    point_list.clean_out();
    tree.nodes_near( point_array[read], point_list );
    while (point_list.size())
    {
      CubitVector* v = point_list.pop();
      int index = v - point_array;
      assert( ((unsigned)index) < pointList.size() );
      if ( (index < read) && ((point_array[read]-*v).length_squared() < tolsqr) )
      {
        index_map[read] = index_map[index];
        break;
      }
    }
    
    if ( index_map[read] == write )
    {
      pointList[write++] = point_array[read];
    }
  }
  pointList.resize(write);
  delete [] point_array;
  
    // Use index_map to update facet list
  IntegerList::iterator itor = triangleData.begin();
  IntegerList::iterator end  = triangleData.end();
  while (itor++ != end)
    for (int count = 3; count--; itor++)
      *itor = index_map[*itor];
  
  delete [] index_map;
  pointsConsolidated = true;
}

void CompSurfFacets::graphics( double tolerance, GMem& gmem )
{
  if (!pointsConsolidated)
  {
    consolidate_points(tolerance);
    pointsConsolidated = true;
  }
  
  gmem.allocate_tri(triangleData.size() / 4);

  assert(triangleData.size() % 4 == 0);
  IntegerList::iterator i_itor = triangleData.begin();
  IntegerList::iterator i_end = triangleData.end();
  int* f_itor = gmem.facet_list();
  
  while (i_itor != i_end)
  {
    i_itor++;
    *(f_itor++) = 3;
    *(f_itor++) = *(i_itor++);
    *(f_itor++) = *(i_itor++);
    *(f_itor++) = *(i_itor++);
  }
  
  PointList::iterator p_itor = pointList.begin();
  PointList::iterator p_end = pointList.end();
  GPoint* g_itor = gmem.point_list();
  
  while (p_itor != p_end)
  {
    g_itor->x = (float)p_itor->x();
    g_itor->y = (float)p_itor->y();
    g_itor->z = (float)p_itor->z();
    g_itor++;
    p_itor++;
  }
  
  gmem.fListCount = triangleData.size();
  gmem.pointListCount = pointList.size();
}

