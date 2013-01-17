/**
 * \file mergechk.cpp
 *
 * \brief mergechk, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * performs imprints between all the bodies, merges them, and writes information
 * on the results.  It also performs pairwise intersections between the
 * bodies to check for overlaps.  Results are written to stardard output.
 *
 */

#undef NDEBUG
#include <cassert>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "InitCGMA.hpp"

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **);
CubitStatus evaluate_overlaps();
CubitStatus imprint_bodies();
CubitStatus print_unmerged_surfaces();
CubitStatus hollow();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
    // Start up CGM
  CubitStatus status = InitCGMA::initialize_cgma("ACIS");
  if (CUBIT_SUCCESS != status) return 1;

  //Do hollow operation to make thick body.
  status = hollow();
  if (status == CUBIT_FAILURE) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val;
  
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, const char **argv) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
    // For each file, open and read the geometry
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    status = gti->import_solid_model(argv[i], "ACIS_SAT");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", argv[i]);
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus hollow()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/hollow.sat";
  CubitStatus status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);
  else if (gti->num_bodies() == 0) {
    PRINT_WARNING("No bodies read; exiting.\n");
    int ret_val = ( CubitMessage::instance()->error_count() );

    exit(ret_val);
  }
  
  //test making thick body.
  DLIList<Body*> new_bodies;
  gti->bodies(new_bodies);
  double d = new_bodies.get()->measure(); //d = 518.3627
  int n = new_bodies.get()->num_ref_faces(); //n = 5
  assert(n == 5);
  assert(d-518.3627 < 0.0001);

  //find the top most surface as the opening of the thick body.
  DLIList<RefFace*> ref_faces;
  new_bodies.get()->ref_faces(ref_faces);
  CubitVector center(0,0,10);
  for(int i = 0; i < n; i++)
  {
    if(ref_faces.step_and_get()->is_planar() &&
       ref_faces.get()->center_point() == center )
      break;
  }
  DLIList<RefFace*> faces_to_remove;
  faces_to_remove.append(ref_faces.get());
  DLIList<Body*> from_bodies;
  from_bodies = new_bodies;
  new_bodies.clean_out();
  CubitStatus stat = gmti->hollow(from_bodies, faces_to_remove, new_bodies, -.2);
  //Created volume(s): 2
  //Destroyed volume(s): 1
  n = new_bodies.get()->num_ref_faces(); //n = 9
  d = new_bodies.get()->measure(); //d = 72.4074
  assert(n == 9);
  assert(d - 72.4074 < 0.0001);
  return stat ;
}

