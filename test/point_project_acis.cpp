/**
 * \file point.project.cpp
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs checks for point projects to trimmed surface. If a point is
 * is projected to empty space, project it to one of the curves of the surface
 * instead.
 */

#undef NDEBUG
#include <cassert>

#include "GeometryModifyEngine.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitMessage.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"
#include "Surface.hpp"
#include "BodySM.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, const char **, bool local = false);
CubitStatus point_project();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");

// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma("ACIS");
  if (CUBIT_SUCCESS != status) return 1;

  //Do make point.
  status = point_project();
  if (status == CUBIT_FAILURE)
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 2 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
    return ret_val;
  }
  return 0;

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
    status = CubitCompat_import_solid_model(filename.c_str(), "ACIS_SAT");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename.c_str());
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus point_project()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();

  const char *argv = "holysurf.sat";
  CubitStatus status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);
  
  DLIList<RefFace*> ref_surfaces;
  gti->ref_faces (ref_surfaces);
  Surface* surf;
  surf = ref_surfaces.get_and_step()->get_surface_ptr();
   
  CubitVector point(0.0,3.,-.65); 
  CubitVector on_surf;
  surf->closest_point_trimmed(point, on_surf); 
  assert (on_surf.z() == 0.5);
  assert (on_surf.y() == 2);
  assert (on_surf.x() == 0 );

  CubitVector p1(0.1, 0.45, -.6 );
  surf->closest_point_trimmed(p1, on_surf);
  assert (on_surf.z() == 0.5);
  assert (on_surf.y() < 0.3882 && on_surf.y() > 0.3881);
  assert (on_surf.x() < 0.2237 && on_surf.x() > 0.2236 );
  return CUBIT_SUCCESS;
}
