

#include "stdio.h"

#include "GMem.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "CubitDefines.h"
#include "CubitBox.hpp"

#include "TestUtilities.hpp"

int GraphicsData(int argc, char* argv[])
{
  GeometryModifyTool::instance()->brick(10,10,10);

  GeometryQueryEngine* engine = GeometryQueryTool::instance()->get_gqe();

  if(!engine)
  {
    printf("no engine\n");
    return 1;
  }
 
  CubitBox uninit_box;
  
  CubitBox model_bbox(
      CubitVector(-5, -5, -5),
      CubitVector(5, 5, 5)
      );


  DLIList<RefEdge*> edges;
  GeometryQueryTool::instance()->ref_edges(edges);

  // expect 24 edges
  if(edges.size() != 12)
  {
    printf("wrong edge count\n");
    return 1;
  }

  // get graphics data for the curves

  CubitBox g_curves = uninit_box;
  GMem gmem;
  for(int i=0; i<edges.size(); i++)
  {
    RefEdge* edge = edges.get_and_step();
    if(engine->get_graphics(edge->get_curve_ptr(), &gmem) != CUBIT_SUCCESS || 
        gmem.pointListCount == 0)
    {
      printf("got no edge facet data for edge %i\n", edge->id());
    }
    else
    {
      for(int j=0; j<gmem.pointListCount; j++)
      {
        const GPoint& pt = gmem.point_list()[j];
        g_curves |= CubitVector(pt.x, pt.y, pt.z);
      }
    }
  }

  printf("curve bbox range = %f %f %f\n",
      g_curves.x_range(), g_curves.y_range(), g_curves.z_range());

  if(!cubit_box_identical(g_curves, model_bbox, GEOMETRY_RESABS))
  {
    printf("didn't get accurate facet data for curves\n");
    return 1;
  }

  
  CubitBox g_surfaces = uninit_box;
  DLIList<RefFace*> faces;
  GeometryQueryTool::instance()->ref_faces(faces);

  // expect 6 edges
  if(faces.size() != 6)
  {
    printf("wrong face count\n");
    return 1;
  }

  // get graphics data for the surfaces
  for(int i=0; i<faces.size(); i++)
  {
    RefFace* face = faces.get_and_step();
    if(engine->get_graphics(face->get_surface_ptr(), &gmem) != CUBIT_SUCCESS || 
        gmem.pointListCount == 0 || gmem.fListCount == 0)
    {
      printf("got no face facet data for face %i\n", face->id());
    }
    else
    {
      const GPoint* pts = gmem.point_list();
      const int* facets = gmem.facet_list();
      int facet_size = gmem.fListCount;
      for(int j=0; j<facet_size;)
      {
        int num = facets[j++]; 
        for(int k=0; k<num; k++)
        {
          const GPoint& pt = pts[facets[j++]];
          g_surfaces |= CubitVector(pt.x, pt.y, pt.z);
        }
      }
    }
  }
  
  printf("surface bbox range = %f %f %f\n",
      g_surfaces.x_range(), g_surfaces.y_range(), g_surfaces.z_range());
  
  if(!cubit_box_identical(g_surfaces, model_bbox, GEOMETRY_RESABS))
  {
    printf("didn't get accurate facet data for surfaces\n");
    return 1;
  }


  return 0;
}

