/**
 * \file imprint_bug.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "Shell.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CastTo.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "OCCShell.hpp"
#include "TopoDS_Shape.hxx"
#include "InitCGMA.cpp"
#include "CubitCompat.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, const char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma("OCC");
  if (CUBIT_SUCCESS != status) return 1;

  //Do make point.
  status = make_Point();
  if (status == CUBIT_FAILURE) 
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

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, const char **argv, bool local)
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    std::string filename( local ? "./" : SRCPATH );
    filename += argv[i];
    status = CubitCompat_import_solid_model(filename.c_str(), "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename.c_str());
      abort();
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus make_Point()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  OCCQueryEngine::instance();
  OCCModifyEngine* ome = OCCModifyEngine::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = "Solid_2.brep";
  CubitStatus  status = read_geometry(1, &argv, true);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 1 volume.

  argv = "Solid_7.brep";
  status = read_geometry(1, &argv, true);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 1 volume.

  DLIList<Body*> from_bodies;
  DLIList<Body*> new_bodies;
  gti->bodies(from_bodies);
  //Bug 247 test, when this bug get fixed, the comment should be removed.
//  status = gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE);

  DLIList<Body*> bodies;
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  DLIList<RefEntity*>  free_entities;
  gti->get_free_ref_entities(free_entities);
  for(int i = 0; i < free_entities.size(); i++)
    gti->delete_RefEntity(free_entities.get_and_step()); 

  Body* body = gmti->brick(10, 10, 10);
  DLIList<Shell*> shells;
  body->shells(shells);
  ShellSM* shell = shells.get()->get_shell_sm_ptr(); 
  OCCShell* occ_shell = CAST_TO(shell, OCCShell);
  DLIList<TopologyBridge*> parents;
  occ_shell->get_parents_virt(parents);
  assert (parents.size() == 1);
  return status;
}
