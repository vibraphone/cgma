/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */
#include "config.h"
#include "CpuTimer.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitUtil.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitObserver.hpp"
#include "CastTo.hpp"
#include "OCCModifyEngine.hpp"
#include "AppUtil.hpp"
#include "RefEntityFactory.hpp"
#include "RefEdge.hpp"
#include "BodySM.hpp"
#include "Lump.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "OCCShell.hpp"
#include "TopoDS_Shape.hxx"
#include "RefEntityName.hpp"
#include "RefEntityFactory.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{

  CubitObserver::init_static_observers();
    // Initialize the GeometryTool
  
  CGMApp::instance()->startup( argc, argv );
  OCCQueryEngine::instance();
  OCCModifyEngine::instance();

    // If there aren't any file arguments, print usage and exit
  //if (argc == 1) {
  //  PRINT_INFO("Usage: mergechk <geom_file> [<geom_file> ...]\n");
  //  exit(0);
  //}
  
  CubitStatus status = CUBIT_SUCCESS;



  //Do make point.
  status = make_Point();
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
    int length = filename.length();
    std::string sub1 = filename.substr(length - 3);
    std::string sub2 = filename.substr(length - 4);
    char * cstr1;
    cstr1 = new char [length+1];
    char * cstr2;
    cstr2 = new char [length+1];
    strcpy (cstr1, sub1.c_str());
    strcpy (cstr2, sub2.c_str());

    if(!strcmp("OCC", cstr1) || 
       !strcmp("occ", cstr1 )||
       !strcmp("brep", cstr2)||
       !strcmp("BREP", cstr2) )
      status = gti->import_solid_model(filename.c_str(), "OCC");
    else if(!strcmp("step", cstr2) ||
            !strcmp("STEP", cstr2)||
            !strcmp("stp", cstr1)  ||
            !strcmp("STP", cstr1))
      status = gti->import_solid_model(filename.c_str(), "STEP");
    else if(!strcmp("iges", cstr2) ||
            !strcmp("IGES", cstr2) ||
            !strcmp("igs", cstr1)  ||
            !strcmp("IGS", cstr1))
      status = gti->import_solid_model(filename.c_str(), "IGES");
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

  DLIList<Body*> bodies;
  DLIList<RefEntity*>  free_entities;

  //Read in the geometry from iges file
  const char *argiges = "ex3.iges";
  CubitStatus status = read_geometry(1, &argiges, false);
  if (status == CUBIT_FAILURE) exit(1);

  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "ex3.occ";
  const char * filetype = "OCC";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  const char *argiges2 = "ex2.iges";
  status = read_geometry(1, &argiges2, false);
  if (status == CUBIT_FAILURE) exit(1);

  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "ex2.occ";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  const char *argstep = "proe.stp";
  status = read_geometry(1, &argstep, false);
  if (status == CUBIT_FAILURE) exit(1);
  
  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "proe.occ";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  // Read in the geometry from files specified on the command line
  const char *argv = "stitch.name_occ";
  status = read_geometry(1, &argv, false);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "beforesub.occ";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->bodies(bodies); 

  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  const char *argv1 = "beforesub.occ";
  status = read_geometry(1, &argv1, true);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  //export the newly read-in file
  filename = "beforesub2.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  DLIList<Body*> new_bodies;
  DLIList<Body*> from_bodies;
  from_bodies.append(bodies.get());
  Body* tool_body = bodies.step_and_get();  
  rsl = gmti->subtract(tool_body,from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  //Created volume(s): 5, 6
  //Destroyed volume(s): 3, 4
  double d = new_bodies.step_and_get()->measure();
  CubitVector v = new_bodies.get()->center_point();
  int n = new_bodies.get()->num_ref_faces();
  // n = 6
  //new bodies has 2 bodies, one has a volume = 10 and the other has a 
  //volume = 50; each of them has 6 ref_faces, of which 3 are new and 3 are
  //remaining (unchanged or modified).

  filename = "aftersub.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 
  
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  const char *argv2 = "unite1.occ";
  status = read_geometry(1, &argv2, false);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  from_bodies.clean_out();
  new_bodies.clean_out();
  gti->bodies(from_bodies);
  status = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  assert(status);

  filename = "unite2.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  const char *argv3 = "unite1.occ";
  status = read_geometry(1, &argv3, false);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  //change the order of the two bodies,and unite, see the united name unchanged.
  new_bodies.clean_out();
  bodies.clean_out();
  gti->bodies(bodies);
  from_bodies.clean_out();
  from_bodies.append(bodies.step_and_get());
  from_bodies.append(bodies.step_and_get());

  status = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  assert(status);
  filename = "unite3.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

    // Read in the geometry from files specified on the command line
  const char *argv4 = "unite4.occ";
  status = read_geometry(1, &argv4, false);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  from_bodies.clean_out();
  new_bodies.clean_out();
  gti->bodies(from_bodies);
  status = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  assert(status);

  filename = "unite5.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  status = read_geometry(1, &argv4, false);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  //change the order of the two bodies, and unite, see the name change.
  new_bodies.clean_out();
  bodies.clean_out();
  gti->bodies(bodies);
  from_bodies.clean_out();
  from_bodies.append(bodies.step_and_get());
  from_bodies.append(bodies.step_and_get());

  status = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  assert(status);
  filename = "unite6.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);
  return CUBIT_SUCCESS;
}
