/**
 * \file AngleCalc.cpp
 *
 * \brief AngleCalc, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>

#include "stdio.h"

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
#include "RefFace.hpp"
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
int AngleCalc();

// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus s = InitCGMA::initialize_cgma( ENGINE );
  if (CUBIT_SUCCESS != s) return 1;

  //Do tests.
  int rsl = AngleCalc();
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

int AngleCalc()
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

