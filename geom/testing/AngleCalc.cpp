

#include "stdio.h"

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
#include "RefFace.hpp"
#include "CubitDefines.h"
#include "CubitBox.hpp"

#include "TestUtilities.hpp"

int AngleCalc(int argc, char* argv[])
{
  Body* brick = GeometryModifyTool::instance()->brick(1,2,4);
  if(!brick)
  {
    printf("failed to make brick\n");
    return 1;
  }

  // compute angles at each vertex and be sure they are 90 degrees
  DLIList<RefFace*> faces;
  GeometryQueryTool::instance()->ref_faces(faces);

  bool errors = false;

  for(int i=0; i<faces.size(); i++)
  {
    // get loops
    DLIList<DLIList<CoEdge*>*> loops;
    RefFace* face = faces.get_and_step();
    face->co_edge_loops(loops);

    for(int j=0; j<loops.size(); j++)
    {
      for(int k=0; k<loops[j]->size(); k++)
      {
        CoEdge* edge1 = loops[j]->next(k);
        CoEdge* edge2 = loops[j]->next(k+1);

        double angle = GeometryQueryTool::instance()->geometric_angle(edge1, edge2);
        angle = angle * (180/CUBIT_PI);
        if(fabs(90 - angle) > 0.01)
        {
          RefEdge* redge1 = edge1->get_ref_edge_ptr();
          RefEdge* redge2 = edge2->get_ref_edge_ptr();
          RefVertex* vert = redge1->common_ref_vertex(redge2);

          printf("wrong angle at vertex (%f) %i with surface %i\n", 
              angle, vert->id(), face->id());
          errors = true;
        }

      }
      delete loops[j];
    }
  }

  return errors ? 1 : 0;
}

