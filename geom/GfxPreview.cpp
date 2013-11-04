

#include "GfxPreview.hpp"
#include "GMem.hpp"
#include "CubitVector.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "GeometryQueryEngine.hpp"
#include "Point.hpp"
#include <assert.h>
#include "DrawingToolDefines.h"

GfxPreview* GfxPreview::mInstance = 0;

GfxPreview::GfxPreview()
{
  if(mInstance)
  {
    assert(0);
    // if you want a new instance in place of a previous one,
    // delete the previous one first
  }
  mInstance = this;
}

GfxPreview::~GfxPreview()
{
  mInstance = 0;
}


// This clears preview data ONLY
void GfxPreview::clear()
{
  if(!mInstance) return;
  mInstance->p_clear();
}

// causes the current window to be redrawn
void GfxPreview::flush()
{
  if(!mInstance) return;
  mInstance->p_flush();
}

// Geometry Drawing Functions
    
void GfxPreview::draw_ref_entity(RefEntity* entity, int color, CubitTransformMatrix* mat)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_entity(entity, color, mat);
}

void GfxPreview::draw_ref_vertex(RefVertex* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_vertex(entity, color);
}

void GfxPreview::draw_ref_edge(RefEdge* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_edge(entity, color);
}

void GfxPreview::draw_ref_face(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face(entity, color);
}

void GfxPreview::draw_ref_face_edges(RefFace* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_face_edges(entity, color);
}

void GfxPreview::draw_ref_volume(RefVolume* entity, int color, int mode)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_volume(entity, color, mode);
}

void GfxPreview::draw_ref_volume_edges(RefVolume* entity, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_volume_edges(entity, color);
}

void GfxPreview::draw_ref_body(Body* entity, int color, int mode)
{
  if(!mInstance) return;
  mInstance->p_draw_ref_body(entity, color, mode);
}

// Generic Primitive Drawing Functions

void GfxPreview::draw_box(CubitBox& box, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_box(box, color);
}

// draw point x,y,z of color
void GfxPreview::draw_point(float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(x, y, z, color);
}

void GfxPreview::draw_point(const CubitVector& vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(vector, color);
}

void GfxPreview::draw_point(CubitVector* vector, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_point(*vector, color);
}

// draw line of points {x,y,z}, {x,y,z} of color
void GfxPreview::draw_line(float x1, float y1, float z1, float x2, float y2, float z2, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x1, y1, z1, x2, y2, z2, color);
}

void GfxPreview::draw_line(const CubitVector& x, const CubitVector& y, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_line(x, y, color);
}

// draw a polyline with a number of points of a color
void GfxPreview::draw_polyline(const GPoint* points, int num_points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_polyline(points, num_points, color);
}

// draw polygons given the coordinates and the connectivity.  Note that the
// connectivity (face_list) is given as the number of points in the polygon 
// followed by the point ids for that polygon.  So, the num_face_points is
// the total number of data points in the face_list such that
// num_face_points = number of points + number of polygons
void GfxPreview::draw_polygons(int num_points, const float* xyzs, int num_face_points, 
                               const int* face_list, int color, int mode )
{
  if (!mInstance) 
    return;
  mInstance->p_draw_polygons( num_points, xyzs, num_face_points, face_list, color, mode);
}

// draw a polyline with a number of points of a color
void GfxPreview::draw_polygon(const GPoint* points, int num_points, int color, int border_color, bool filled)
{
  if(!mInstance) return;
  mInstance->p_draw_polygon(points, num_points, color, border_color, filled);
}

// draw a quad of a color
void GfxPreview::draw_quad(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_quad(points, color);
}

// draw a tri of a color
void GfxPreview::draw_tri(const GPoint* points, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_tri(points, color);
}

void GfxPreview::draw_vector(const CubitVector& tail, const CubitVector& head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(tail, head, color);
}

void GfxPreview::draw_vector(const GPoint* tail, const GPoint* head, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_vector(CubitVector(tail->x, tail->y, tail->z), CubitVector(head->x, head->y, head->z), color);
}

// draw a label at position x,y,z of color
void GfxPreview::draw_label(const char* label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}
void GfxPreview::draw_label(int label, float x, float y, float z, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_label(label, x, y, z, color);
}

void GfxPreview:: draw_cylinder(const CubitVector& axis, const CubitVector& origin, 
                                 CubitBox& bounding_box, float radius, int color)
{
  if(!mInstance) return;
  mInstance->p_draw_cylinder(axis, origin, bounding_box, radius, color);
}

void GfxPreview:: draw_frustum( const CubitVector axis , const CubitVector origin1, const CubitVector origin2,double dRad1 , double dRad2,
                               CubitBox& bounding_box ,int color)
{
  if(!mInstance) return;
  mInstance->p_draw_frustum( axis , origin1, origin2, dRad1 , dRad2 , bounding_box, color);
}

void GfxPreview::draw_point(TBPoint *pt, int color)
{
  CubitVector v = pt->coordinates();
  draw_point(v, color);
  flush();
}

void GfxPreview::draw_curve(Curve *curve, int color)
{
  int num_divs = 20;
  double lo, hi;
  curve->get_param_range(lo, hi);
  double dv = (hi-lo)/(double)num_divs;
  int i;
  double param = lo;
  for(i=0; i<num_divs; i++)
  {
    CubitVector p1, p2;
    curve->position_from_u(param, p1);
    param += dv;
    curve->position_from_u(param, p2);
    draw_line(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), color);
  }
  flush();
  DLIList<TBPoint*> pts;
  curve->points(pts);
  for(i=pts.size(); i>0; i--)
  {
    TBPoint *cur_pt = pts.get_and_step();
    draw_point(cur_pt, color);
  }
}

void GfxPreview::draw_surface(Surface *surf, int color)
{
  int num_divs = 20;
  double ulo, uhi, vlo, vhi;
  surf->get_param_range_U(ulo, uhi);
  surf->get_param_range_V(vlo, vhi);
  double du = (uhi-ulo)/(double)num_divs;
  double dv = (vhi-vlo)/(double)num_divs;
  int i, j;
  double uparam, vparam;
  uparam = ulo;
  for(i=0; i<num_divs; i++)
  {
    vparam = vlo;
    for(j=0; j<num_divs; j++)
    {
      CubitVector p1, p2;
      p1 = surf->position_from_u_v(uparam, vparam);
      vparam += dv;
      p2 = surf->position_from_u_v(uparam, vparam);
      draw_line(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), color);
    }
    uparam += du;
  }
  vparam = vlo;
  for(i=0; i<num_divs; i++)
  {
    uparam = ulo;
    for(j=0; j<num_divs; j++)
    {
      CubitVector p1, p2;
      p1 = surf->position_from_u_v(uparam, vparam);
      uparam += du;
      p2 = surf->position_from_u_v(uparam, vparam);
      draw_line(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), color);
    }
    vparam += dv;
  }
  flush();
}

void GfxPreview::draw_surface_facets_shaded(Surface *surf, int color)
{
  GMem g_mem;

  surf->get_geometry_query_engine()->get_graphics(surf, &g_mem);

  const float* xyzs = reinterpret_cast<const float*>(g_mem.point_list());
  GfxPreview::draw_polygons(g_mem.pointListCount, xyzs, 
                            g_mem.fListCount, g_mem.facet_list(),
                            color, CUBIT_SMOOTH_SHADING);
  flush();
}

