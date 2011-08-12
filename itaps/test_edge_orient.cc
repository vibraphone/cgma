#include "iGeom.h"
#include <iostream>
#define CHECK( STR ) if (err != iBase_SUCCESS) return print_error( STR, err, geom, __FILE__, __LINE__ )

#define FILE_NAME1  "ilc_problem_surf8.stp"
#define FILE_NAME2  "brick_2.stp"

static bool print_error( const char* desc, 
                         int err,
                         iGeom_Instance geom,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';
  
  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;
  
  return false; // must always return false or CHECK macro will break
}

int main(int argc, char *argv[])
{
  // initialize the Mesh
  int i, j, err;
  iGeom_Instance geom;
  iGeom_newGeom(NULL, &geom, &err, 0);
  CHECK( "ERROR : can not make a new iGeom instance." );

  // read in the geometry
  std::string input_file;
  input_file += FILE_NAME1;
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry." );

  input_file = FILE_NAME2;
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry." );

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet(geom, &root_set, &err);
  CHECK("Failed to get root set.");
  
  // get a volume
  iBase_EntityHandle* vols = NULL;
  int v_alloc = 0;
  int v_size = 0;
  iGeom_getEntities(geom, root_set, iBase_REGION, &vols,
                    &v_alloc, &v_size, &err);
  CHECK("Failed to get volumes.");

  // get faces of volume
  iBase_EntityHandle* all_faces = NULL;
  int af_alloc = 0;
  int af_size = 0;
  iGeom_getEntAdj(geom, vols[0], iBase_FACE, &all_faces, &af_alloc,
                &af_size, &err);
  CHECK("Failed to get faces.");
  
  // CHECK 1 : all face senses are FORWARD respect to parent volume
  for (i = 0; i < af_size; i++) { // for all faces
    // get face sense compared by volume
    int face_sense;
    iGeom_getEntNrmlSense(geom, all_faces[i], vols[0], &face_sense, &err);
    CHECK("Failed to get face sense.");
    
    // check if face sense is FORWARD respect to parent volume
    if (face_sense != 1) {
      std::cerr << "Error: face sense is not FORWARD." << std::endl;
      return 1;
    }
  }

  // get all edges
  iBase_EntityHandle* edges = NULL;
  int e_alloc = 0;
  int e_size = 0;
  iGeom_getEntAdj(geom, vols[0], iBase_EDGE, &edges, &e_alloc,
                  &e_size, &err);
  CHECK("Failed to get edges.");

  // CHECK 2 : Edge senses to 2 parent faces should be opposite
  for (i = 0; i < e_size; i++) { // for all edges
    // get parent faces 
    iBase_EntityHandle* faces = NULL;
    int f_alloc = 0;
    int f_size = 0;
    iGeom_getEntAdj(geom, edges[i], iBase_FACE, &faces, &f_alloc,
                    &f_size, &err);
    CHECK("Failed to get edges.");

    // check if # of parent faces of edges is 2
    if (f_size != 2) {
      std::cerr << "Error: # of parent faces of brick edges should be 2." << std::endl;
      return 1;
    }

    // get edge senses respect to parent faces
    int sense[2];
    for (j = 0; j < 2; j++) {
      iGeom_getEgFcSense(geom, edges[i], faces[j], &sense[j], &err);
      CHECK("Failed to get edge sense.");
    }
    
    // check if edge senses are opposite
    if (sense[0]*sense[1] != -1) {
      std::cerr << "Error: Edge senses to 2 parent faces should be opposite." << std::endl;
      return 1;
    }
  }

  std::cout << "All tests are passed." << std::endl;

  return 0;
}
