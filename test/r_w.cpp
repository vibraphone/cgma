/**
 * \file makept.cpp
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
#include "CubitMessage.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"

#include <algorithm>

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **, bool local = false,  
                          const char* unit = NULL);
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
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val;
  
}

std::string type_from_file_name( const std::string& filename )
{
  size_t dot_pos = filename.find_last_of( '.' );
  if (dot_pos == std::string::npos)
    return std::string();
 
  std::string extension = filename.substr( dot_pos + 1 );
  std::transform( extension.begin(), extension.end(), extension.begin(), tolower );
  if (extension == "occ" || extension == "brep")
    return "OCC";
  else if (extension == "step" || extension == "stp")
    return "STEP";
  else if (extension == "iges" || extension == "igs")
    return "IGES";
  else
    return std::string();
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, const char **argv, bool local,
                          const char* unit ) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    std::string type = type_from_file_name( argv[i] );
    if (type.empty()) // just guess OCC
      type = "OCC";
    std::string filename( local ? "./" : SRCPATH );
    char const* local_name = argv[i];
    filename += local_name;  
    status = CubitCompat_import_solid_model(filename.c_str(), type.c_str(),
                                            NULL, unit);
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

  DLIList<Body*> bodies;
  DLIList<RefEntity*>  free_entities;

  //Read in the geometry from iges file
  const char *argiges2 = "diffuser.iges";
  CubitStatus status = read_geometry(1, &argiges2, false, "IN");
  //Constructed 7 Volumes: 1 to 7
  if (status == CUBIT_FAILURE) exit(1);

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "diffuser.occ";
  const char * filetype = "OCC";

  CubitStatus rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  //Exported:   7 OCC Bodies to diffuser.occ
  assert(num_ents_exported == 7);

  filetype = "IGES";
  filename = "diffuser_exp.iges";
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version, NULL, "IN");

  //Exported:  7 OCC Curves to diffuer_exp.iges
  assert(num_ents_exported == 7);


  //check that the unit is exported as 1HM
  std::ifstream Myfile;
  std::string line;
  std::string search = "1,4HINCH";
  Myfile.open ("diffuer_exp.iges");
  int found = 0;
  size_t offset;
  if(Myfile.is_open())
  {
    while(!Myfile.eof())
    {
      getline(Myfile,line);
      if ((offset = line.find(search, 0)) != std::string::npos)
        found ++;
    }
    Myfile.close();
  }

  assert (found == 1);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert (free_entities.size() == 0);

  const char *argiges = "ex3.iges";
  status = read_geometry(1, &argiges, false, NULL);
  //Constructed 18 Free Curves: 1 to 18
  if (status == CUBIT_FAILURE) exit(1);

  rsl = CUBIT_SUCCESS;
  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "ex3.occ";
  filetype = "OCC";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  //Exported:  18 OCC Curves to ex3.occ
  assert(num_ents_exported == 18);

  filetype = "IGES";
  filename = "ex3export.iges";
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  //Exported:  18 OCC Curves to ex3export.iges
  assert(num_ents_exported == 18);

  filetype = "STEP";
  filename = "ex3export.step";
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  //Exported:  18 OCC Curves to ex3export.step
  assert(num_ents_exported == 18);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  const char *argiges3 = "cduct.igs";
  status = read_geometry(1, &argiges2, false, "M");
  //Constructed 7 Volumes: 8 to 14
  if (status == CUBIT_FAILURE) exit(1);

  filetype = "IGES";
  filename = "cduct_exp.iges";
  ref_entity_list.clean_out();
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version, NULL, "M");

  //7 OCC Bodies to cduct_exp.iges
  assert(num_ents_exported == 7);


  //check that the unit is exported as 1HM
  search = "6,1HM";
  Myfile.open ("cduct_exp.iges");
  found = 0;
  if(Myfile.is_open())
  {
    while(!Myfile.eof())
    {
      getline(Myfile,line);
      if ((offset = line.find(search, 0)) != std::string::npos)
        found ++;
    }
    Myfile.close();
  }

  assert (found == 1);

  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);

  const char *argstep = "proe.stp";
  status = read_geometry(1, &argstep, false, NULL);
  //Constructed 12 Volumes: 8 to 19
  if (status == CUBIT_FAILURE) exit(1);
  
  filetype = "OCC";
  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "proe.occ";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  //Exported:  12 OCC Bodies to proe.occ
  assert(num_ents_exported == 12);

  filetype = "IGES";
  filename = "proeexport.iges";
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  //Exported:  12 OCC Bodies to proeexport.iges
  assert(num_ents_exported == 12);

  filetype = "STEP";
  filename = "proeexport.step";
  num_ents_exported=0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  //Exported:  12 OCC Bodies to proeexport.step
  assert(num_ents_exported == 12);

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
  status = read_geometry(1, &argv, false, NULL);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  filetype = "OCC";
  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "beforesub.occ";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  assert(num_ents_exported == 2);

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
  status = read_geometry(1, &argv1, true, NULL);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  //export the newly read-in file
  filename = "beforesub2.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  assert(num_ents_exported == 2);
  bodies.clean_out();
  gti->bodies(bodies);
  DLIList<Body*> new_bodies;
  DLIList<Body*> from_bodies;
  from_bodies.append(bodies.get());
  Body* tool_body = bodies.step_and_get();  
  rsl = gmti->subtract(tool_body,from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  //Created volume(s): 24, 25
  //Destroyed volume(s): 22, 23
  double d = new_bodies.step_and_get()->measure();
  CubitVector v = new_bodies.get()->center_point();
  int n = new_bodies.get()->num_ref_faces();
  assert (n == 12);
  CubitVector test_v(5, 5,0.5);
  assert (v == test_v);
  assert (60-d < 0.000000001);
  //new body has 2 volumes, one has a volume = 10 and the other has a 
  //volume = 50; each of them has 6 ref_faces, of which 3 are new and 3 are
  //remaining (unchanged or modified).

  filename = "aftersub.occ";
  ref_entity_list.clean_out();
  num_ents_exported = 0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  assert(num_ents_exported == 1);
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 
  
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  const char *argv2 = "unite1.occ";
  status = read_geometry(1, &argv2, false, NULL);
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
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  assert(num_ents_exported == 1);
  //name attributes in unite2.occ is 
  //CGM_ATTRIB 11 ENTITY_NAME 2 8 volume B 8 volume_B 0 0
  //so it changed to volume_B, and code needs use the second name attrib.
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  const char *argv3 = "unite1.occ";
  status = read_geometry(1, &argv3, false, NULL);
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
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 1);
  //CGM_ATTRIB 11 ENTITY_NAME 2 8 volume A 8 volume_A 0 0 in unite3.occ

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

    // Read in the geometry from files specified on the command line
  const char *argv4 = "unite4.occ";
  status = read_geometry(1, &argv4, false, NULL);
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
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  assert(num_ents_exported == 1);
  //CGM_ATTRIB 11 ENTITY_NAME 1 7 volumeA 0 0 in unite5.occ
 
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  // Read in the geometry from files specified on the command line
  status = read_geometry(1, &argv4, false, NULL);
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
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 1);
 
  //CGM_ATTRIB 11 ENTITY_NAME 1 7 volumeB 0 0 in unite6.occ

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  //testing for copy command which will copy attributes.
  //const char *argv2 = "unite1.occ";
  status = read_geometry(1, &argv2, false, NULL);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  from_bodies.clean_out();
  gti->bodies(from_bodies);
  Body* orig_body = from_bodies.get();
  Body* copy_body = gmti->copy_body(orig_body);
  assert (copy_body);
  filename = "copy_attrib.occ";
  ref_entity_list.clean_out();
  ref_entity_list.append(copy_body);
  num_ents_exported = 0;
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 1);
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);

  return CUBIT_SUCCESS;
}
