/**
 * \file CreateGeometry.cpp
 *
 * \brief CreateGeometry, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>

#include "stdio.h"

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "Body.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitDefines.h"
#include "CubitBox.hpp"
#include "InitCGMA.hpp"
#include "TestUtilities.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#ifdef TEST_ACIS
#  define ENGINE "ACIS"
#elif defined (TEST_OCC)
#  define ENGINE "OCC"
#else
#  error "Which engine to test?"
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
int CreateGeometry();

// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus s = InitCGMA::initialize_cgma( ENGINE );
  if (CUBIT_SUCCESS != s) return 1;

  //Do tests.
  int rsl = CreateGeometry();
  if (rsl == 1) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val != 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  else
    ret_val = 0;

  return ret_val;
  
}

static int test_brick()
{
  Body* brick = GeometryModifyTool::instance()->brick(1,2,4);
  if(!brick)
  {
    printf("failed to make brick\n");
    return 1;
  }
  if(!cubit_box_identical(brick->bounding_box(),
        CubitBox(CubitVector(-0.5,-1,-2), CubitVector(0.5,1,2)),
        GEOMETRY_RESABS*2.0, true))
  {
    printf("boxes not identical\n");
    return 1;
  }
  return 0;
}

static int test_oriented_brick()
{
  CubitVector center(1, 2, 4);
  CubitVector axes[3] = {
    CubitVector(0.707, 0.707, 0),
    CubitVector(-0.707, 0.707, 0),
    CubitVector(0, 0, 1)
  };
  CubitVector extension(3, 5, 7);

  Body* brick = GeometryModifyTool::instance()->brick(center, axes, extension);
  if(!brick)
  {
    printf("failed to make brick\n");
    return 1;
  }
  CubitBox comp_box(CubitVector(-4.656854, -3.656854, -3.),
          CubitVector(6.656854, 7.656854, 11.));
  CubitBox bnd_box = brick->bounding_box();

  bool identical =  cubit_box_identical(bnd_box, comp_box, GEOMETRY_RESABS*2.0, true);
  if (identical)
    return 0;

  if( bnd_box < comp_box || bnd_box > comp_box*1.09)
  {
    printf("boxes not identical\n");
    return 1;
  }
  return 0;
}

static int test_sphere()
{
  Body* sphere = GeometryModifyTool::instance()->sphere(1);
  if(!sphere)
  {
    printf("failed to make sphere\n");
    return 1;
  }
  CubitBox comp_box( CubitVector(-1,-1,-1), CubitVector(1,1,1));
  CubitBox bnd_box = sphere->bounding_box();

  bool identical =  cubit_box_identical(bnd_box, comp_box, GEOMETRY_RESABS*2.0, true);
  if (identical)
    return 0;

  if( bnd_box < comp_box || bnd_box > comp_box*1.09)
  {
    printf("boxes not identical\n");
    return 1;
  }
  return 0;
}

static int test_torus()
{
  printf("making torus\n");
  Body* torus = GeometryModifyTool::instance()->torus(1, .2);
  if(!torus)
  {
    printf("failed to make torus\n");
    return 1;
  }
  CubitBox comp_box(CubitVector(-1.2,-1.2,-0.2), CubitVector(1.2,1.2,0.2));
  CubitBox bnd_box = torus->bounding_box();

  bool identical =  cubit_box_identical(bnd_box, comp_box, GEOMETRY_RESABS*2.0, true);
  if (identical)
    return 0;

  if( bnd_box < comp_box || bnd_box > comp_box*1.09)
  {
    printf("boxes not identical\n");
    return 1;
  }
  return 0;
}

static int test_planar_sheet()
{
  CubitVector axes[2] = {
    CubitVector(1,0,0),
    CubitVector(.1,.9,0)
    };

  Body* body = GeometryModifyTool::instance()->planar_sheet(CubitVector(1,1,1),
      axes, 2, 3);

  if(!body)
  {
    printf("failed to make planar sheet\n");
    return 1;
  }

  CubitBox comp_box(CubitVector(-0.165647,-0.490826,1.0),
        CubitVector(2.165647,2.490826,1.0));
  CubitBox bnd_box = body->bounding_box();

  bool identical =  cubit_box_identical(bnd_box, comp_box, GEOMETRY_RESABS*2.0, true);
  if (identical)
    return 0;

  if( bnd_box < comp_box || bnd_box > comp_box*1.09)
  {
    printf("boxes not identical\n");
    return 1;
  }


  body = GeometryModifyTool::instance()->planar_sheet(
      CubitVector(0,0,0),
      CubitVector(1,0,0),
      CubitVector(1,1,0),
      CubitVector(0,1,.2)
      );

  if(body)
  {
    printf("should have failed to make planar sheet out of non-planar points\n");
    return 1;
  }

  body = GeometryModifyTool::instance()->planar_sheet(
      CubitVector(0,0,0),
      CubitVector(0,0,0),
      CubitVector(1,1,0),
      CubitVector(0,1,0)
      );

  if(body)
  {
    printf("should have failed to make planar sheet with coincident input points\n");
    return 1;
  }

  return 0;
}

static int test_arc()
{
  RefVertex* pt1 = GeometryModifyTool::instance()->make_RefVertex(CubitVector(0,0,0));
  RefVertex* pt2 = GeometryModifyTool::instance()->make_RefVertex(CubitVector(1,1,0));
  RefVertex* pt3 = GeometryModifyTool::instance()->make_RefVertex(CubitVector(2,0,0));

  RefEdge* arc = GeometryModifyTool::instance()->create_arc_three(pt1, pt2, pt3, false);
  if(!arc)
  {
    printf("failed to make arc\n");
    return 1;
  }

  CubitBox comp_box(CubitVector(0, 0, 0), CubitVector(2, 1, 0));
  CubitBox bnd_box = arc->bounding_box();

  bool identical =  cubit_box_identical(bnd_box, comp_box, GEOMETRY_RESABS*2.0, true);
  if (!identical)
  {
    if( bnd_box < comp_box || bnd_box > comp_box*1.09)
    {
      printf("boxes not identical\n");
      return 1;
    }
  }

  // previously free curves at end points must be consumed
  DLIList<RefVertex*> all_verts;
  GeometryQueryTool::instance()->ref_vertices(all_verts);
  if(all_verts.size() != 3)
  {
    printf("vertices not consumed properly with curve creation\n");
    return 1;
  }

  GeometryQueryTool::instance()->delete_RefVertex(pt2);

  RefVertex* pt4 = GeometryModifyTool::instance()->make_RefVertex(CubitVector(1,-1,0));
  RefEdge* arc2 = GeometryModifyTool::instance()->create_arc_three(pt1, pt4, pt3, true);

  CubitBox comp_box2(CubitVector(0, -1, 0), CubitVector(2, 1, 0));
  bnd_box = arc2->bounding_box();

  identical =  cubit_box_identical(bnd_box, comp_box2, GEOMETRY_RESABS*2.0, true);
  if (identical)
    return 0;

  if( bnd_box < comp_box2 || bnd_box > comp_box2*1.09)
  {
    printf("boxes not identical\n");
    return 1;
  }

  all_verts.clean_out();
  GeometryQueryTool::instance()->ref_vertices(all_verts);
  if(all_verts.size() != 4)
  {
    printf("vertices not consumed properly with curve creation\n");
    return 1;
  }

  return 0;
}

int CreateGeometry()
{
  int ret = 0;

  ret |= test_brick();
  GeometryQueryTool::instance()->delete_geometry();

  ret |= test_oriented_brick();
  GeometryQueryTool::instance()->delete_geometry();

  ret |= test_sphere();
  GeometryQueryTool::instance()->delete_geometry();

  ret |= test_torus();
  GeometryQueryTool::instance()->delete_geometry();

  ret |= test_planar_sheet();
  GeometryQueryTool::instance()->delete_geometry();

  ret |= test_arc();
  GeometryQueryTool::instance()->delete_geometry();

  return ret;
}

