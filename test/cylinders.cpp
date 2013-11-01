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
#include <string>
#include <cctype>

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

#ifdef HAVE_OCC
#define FILENAME "s5.stp" 
#else
#define FILENAME  "s5.sat"
#endif
// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");

int findString(const char *filename, std::string search)
{
  std::ifstream Myfile;
  Myfile.open (filename);
  int found = 0;
  std::string line;
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
  return found;
}
// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma();
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
  std::transform( extension.begin(), extension.end(), extension.begin(), ::tolower );
  if (extension == "occ" || extension == "brep")
    return "OCC";
  else if (extension == "step" || extension == "stp")
    return "STEP";
  else if (extension == "iges" || extension == "igs")
    return "IGES";
  else if (extension == "sat")
    return "ACIS_SAT";
  else
    return std::string();
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
    std::string type = type_from_file_name( argv[i] );
    if (type.empty()) // just guess OCC
      type = "OCC";
    std::string filename( local ? "./" : SRCPATH );
    char const* local_name = argv[i];
    filename += local_name;  
    status = CubitCompat_import_solid_model(filename.c_str(), type.c_str());
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

  DLIList<Body*> bodies, single_body, all_bodies, neighbor_list, new_bodies;
  DLIList<RefEntity*>  free_entities;

  const char *argstep = FILENAME;
  CubitStatus stat = read_geometry(1, &argstep, false);
  //Constructed 12 Volumes: 8 to 19
  if (stat == CUBIT_FAILURE) exit(1);
  
  //get all the bodies, use bodies 1-6 to unite and use united body to 
  //subtract from body 7. 
  gti->bodies(bodies);
  Body* brick = bodies.pop();
  all_bodies.append(brick);

  stat = gmti->unite(bodies, single_body, CUBIT_FALSE);
  assert(CUBIT_SUCCESS == stat && single_body.size() == 1);

  stat = gmti->subtract(single_body.get(), all_bodies, new_bodies, CUBIT_FALSE,CUBIT_FALSE);
  assert(CUBIT_SUCCESS == stat);
 
  std::cout << "Number of resulting bodies = " << new_bodies.size() << std::endl;
  assert(new_bodies.size() == 1);
  //delete all entities
  bodies.clean_out();
  gti->bodies(bodies);
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);
  return CUBIT_SUCCESS;
}
