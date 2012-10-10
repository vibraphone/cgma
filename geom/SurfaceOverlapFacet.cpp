//- Class: SurfaceOverlapFacet
//- Description: Facet definition class for efficient processing
//-              for SurfaceOverlapTool.  
//- Owner: Steve Storm
//- Created: January 26, 2003

#include "SurfaceOverlapFacet.hpp"
#include "GfxPreview.hpp"

AnalyticGeometryTool* SurfaceOverlapFacet::agt = AnalyticGeometryTool::instance();

// Constructor
SurfaceOverlapFacet::SurfaceOverlapFacet( GPoint p[3] )
{
  t.b.x = p[0].x;
  t.b.y = p[0].y;
  t.b.z = p[0].z;

  t.e0.x = p[1].x - p[0].x;
  t.e0.y = p[1].y - p[0].y;
  t.e0.z = p[1].z - p[0].z;

  t.e1.x = p[2].x - p[0].x;
  t.e1.y = p[2].y - p[0].y;
  t.e1.z = p[2].z - p[0].z;

  CubitVector min;
   
  min.set( CUBIT_MIN_3( p[0].x, p[1].x, p[2].x ),
    CUBIT_MIN_3( p[0].y, p[1].y, p[2].y ),
    CUBIT_MIN_3( p[0].z, p[1].z, p[2].z ) );
  
  CubitVector max;
  
  max.set( CUBIT_MAX_3( p[0].x, p[1].x, p[2].x ),
    CUBIT_MAX_3( p[0].y, p[1].y, p[2].y ),
    CUBIT_MAX_3( p[0].z, p[1].z, p[2].z ) );

  boundingBox.reset( min, max );
}

// Destructor
SurfaceOverlapFacet::~SurfaceOverlapFacet()
{
}

double SurfaceOverlapFacet::distance( SurfaceOverlapFacet &other_facet )
{
  double s, t2, u, v;
  return agt->MinTriangleTriangle( t, other_facet.t, s, t2, u, v );
}

double SurfaceOverlapFacet::perimeter()
{
  CubitVector pt0( t.b.x, t.b.y, t.b.z );
  CubitVector pt1( t.b.x + t.e0.x,
                   t.b.y + t.e0.y,
                   t.b.z + t.e0.z );
  CubitVector pt2( t.b.x + t.e1.x,
                   t.b.y + t.e1.y,
                   t.b.z + t.e1.z ); 

  double total_dist = pt0.distance_between(pt1);
  total_dist += pt1.distance_between( pt2 );
  total_dist += pt2.distance_between( pt0 );

  return total_dist;
}

bool SurfaceOverlapFacet::facet_points_within_tol( SurfaceOverlapFacet *other_face, double tolerance )
{
  CubitVector tmp_pt( t.b.x, t.b.y, t.b.z );
  if( this->distance_from_position( tmp_pt ) > tolerance )
    return false;

  tmp_pt.set( t.b.x + t.e0.x,
              t.b.y + t.e0.y,
              t.b.z + t.e0.z );

  if( this->distance_from_position( tmp_pt ) > tolerance )
    return false;

  tmp_pt.set( t.b.x + t.e1.x,
              t.b.y + t.e1.y,
              t.b.z + t.e1.z );

  if( this->distance_from_position( tmp_pt ) > tolerance )
    return false;

  return true;
}

double SurfaceOverlapFacet::distance_from_position( CubitVector &position )
{
  double s,t;  
  Point3 tmp_point;
  tmp_point.x = position.x();
  tmp_point.y = position.y();
  tmp_point.z = position.z();
  return agt->MinPointTriangle( tmp_point, this->t, s, t );
} 

CubitBoolean
SurfaceOverlapFacet::facing( SurfaceOverlapFacet &other_facet )
{
    double norm1[3];
    agt->Normal( t, norm1 );
    double norm2[3];
    agt->Normal( other_facet.t, norm2 );

    // move to the origin
    double o1[3];
    o1[0] = other_facet.t.b.x-t.b.x;
    o1[1] = other_facet.t.b.y-t.b.y;
    o1[2] = other_facet.t.b.z-t.b.z;

    double dot_p = agt->dot_vec(norm1,o1);
    return (CubitBoolean)(dot_p >= 0.0);
}

double
SurfaceOverlapFacet::angle( SurfaceOverlapFacet &other_facet )
{
  return agt->Angle( t, other_facet.t ) * RADtoDEG;
}

double
SurfaceOverlapFacet::projected_overlap( SurfaceOverlapFacet &other_facet, CubitBoolean draw_overlap )
{
  double tmp_double = agt->ProjectedOverlap( t, other_facet.t, draw_overlap );

  if( tmp_double > 0.00 ) 
  {
    CubitVector edge0(t.e0.x, t.e0.y, t.e0.z);  
    CubitVector edge1(t.e1.x, t.e1.y, t.e1.z);  
    CubitVector normal = edge0 * edge1;
    double area_facet1 = normal.length() / 2;

    edge0.set(other_facet.t.e0.x, other_facet.t.e0.y, other_facet.t.e0.z);  
    edge1.set(other_facet.t.e1.x, other_facet.t.e1.y, other_facet.t.e1.z);  
    normal = edge0 * edge1;
    double area_facet2 = normal.length() / 2;
    
    //don't report overlapping area between facets unless it is greater 
    //than one hundredth of the area of the smaller facet
    if( area_facet1 < area_facet2 )
    {
      if( tmp_double < (area_facet1*0.01))
        tmp_double = 0.0;
    }
    else if( tmp_double < (area_facet2*0.01 ))
      tmp_double = 0.0;
  }
  return tmp_double;
}


void SurfaceOverlapFacet::draw( int color ) 
{
  CubitVector point1 = CubitVector( t.b.x, t.b.y, t.b.z );    
  CubitVector point2 = CubitVector( t.e0.x + t.b.x, 
                                    t.e0.y + t.b.y, 
                                    t.e0.z + t.b.z );    
  CubitVector point3 = CubitVector( t.e1.x + t.b.x, 
                                    t.e1.y + t.b.y, 
                                    t.e1.z + t.b.z );    

  GfxPreview::draw_line( point1, point2, color );
  GfxPreview::draw_line( point2, point3, color );
  GfxPreview::draw_line( point1, point3, color );
  
  return;
}

CubitVector SurfaceOverlapFacet::centroid()
{

  CubitVector point1 = CubitVector( t.b.x, t.b.y, t.b.z );    
  CubitVector point2 = CubitVector( t.e0.x + t.b.x, 
                                    t.e0.y + t.b.y, 
                                    t.e0.z + t.b.z );    
  CubitVector point3 = CubitVector( t.e1.x + t.b.x, 
                                    t.e1.y + t.b.y, 
                                    t.e1.z + t.b.z );    

  point1 += point2;
  point1 += point3;
  point1 /= 3;
  
  return point1;
}

CubitVector SurfaceOverlapFacet::smallest_edge_midpoint()
{
  CubitVector point1 = CubitVector( t.b.x, t.b.y, t.b.z );    
  CubitVector point2 = CubitVector( t.e0.x + t.b.x, 
                                    t.e0.y + t.b.y, 
                                    t.e0.z + t.b.z );    
  CubitVector point3 = CubitVector( t.e1.x + t.b.x, 
                                    t.e1.y + t.b.y, 
                                    t.e1.z + t.b.z );    
  
  double len_12 = point1.distance_between( point2 );
  double len_23 = point2.distance_between( point3 );
  double len_13 = point1.distance_between( point3 );

  if( (len_12 < len_23) && (len_12 < len_13) )
  {
    point1 += point2;
    point1 /= 2;
    return point1;
  }
  else if( len_23 < len_13 )
  {
    point2 += point3;
    point2 /= 2;
    return point2;
  }
  else
  {
    point1 += point3;
    point1 /= 2;
    return point1;
  }
}
