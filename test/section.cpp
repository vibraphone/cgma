/**
 * \file section.cpp
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

  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  OCCModifyEngine* ome = OCCModifyEngine::instance();

  //test for tweak fillet and chamfer
  Body* body = gmti->sphere(53.726642);
  Body* body2 = gmti->brick(186, 372, 372);
  CubitVector move(96, 0, 1);
  gti->translate(body2, move);
  CubitBox box1;
  box1 = body2->bounding_box();
  
  DLIList<Body*> from_bodies, new_bodies;
  from_bodies.append(body);
  CubitStatus rsl = gmti->subtract(body2, from_bodies, new_bodies,
                       CUBIT_FALSE, CUBIT_FALSE);

  int i = new_bodies.size();

  int num_ents_exported=0;
  DLIList<RefEntity*> ref_entities;
  const CubitString cubit_version="10.2";
  const char * filename = "section.occ";
  const char * filetype = "OCC";

  rsl = CubitCompat_export_solid_model(ref_entities, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 1);

  DLIList<Body*> bodies;
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  remove(filename);
  return CUBIT_SUCCESS;
}
