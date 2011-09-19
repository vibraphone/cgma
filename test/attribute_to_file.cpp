#include <iostream>
#include <fstream>
#include <string>

#include "InitCGMA.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "CGMApp.hpp"
#include "CubitAttribManager.hpp"
#include "CAEntityId.hpp"

#include <stdio.h>
#include <stdlib.h>

#ifdef TEST_ACIS
#  define ENGINE "ACIS"
#  define FORMAT "ACIS_SAT"
#elif defined (TEST_OCC)
#  define ENGINE "OCC"
#  define FORMAT "OCC"
#elif defined (TEST_FACET)
#  define ENGINE "FACET"
#  define FORMAT "FACET"
#else
#  error "Which engine to test?"
#endif

#define ASSERT(A) if (!(A)) failed(#A,__FILE__,__LINE__)
CubitBoolean is_files_same(std::string filename1,
                           std::string filename2);

void failed( const char* A, const char* FILE, int LINE )
{
  printf( "Condition failed at %s:%d : %s\n", FILE, LINE, A );
  abort();
}

int main()
{
  CubitStatus s = InitCGMA::initialize_cgma( ENGINE );
  ASSERT(s);
  
  // actuate CA_BODIES attribute and turn on auto flag
  CGMApp::instance()->attrib_manager()->register_attrib_type(CA_ENTITY_ID, "id", "ENTITY_ID",
							     CAEntityId_creator, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->auto_flag(CUBIT_TRUE);

  // make 2 bricks
  int n = 2;
  Body** bricks = new Body*[n];
  DLIList<RefEntity*> export_list;
  for (int i = 0; i < n; i++) {
    bricks[i] = GeometryModifyTool::instance()->brick( 1, 1, 1 );
    ASSERT(bricks[i]);
    s = GeometryQueryTool::instance()->translate( bricks[i], CubitVector(i,0,0) );
    ASSERT(s);
    export_list.append( bricks[i] );
  }

  // export as file
  int junk;
  s = GeometryQueryTool::instance()->export_solid_model( export_list, "bricks2.occ",
                                                         FORMAT, junk, CubitString(__FILE__) );
  ASSERT(s);

  //check that the two single volume bodys' attributes are exported as SINGLELUMP%
  std::ifstream Myfile;
  std::string line;
  std::string search = "SINGLELUMP%";
  Myfile.open ("bricks2.occ");
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
   
  assert (found == 2); 

  s = GeometryQueryTool::instance()->export_solid_model( export_list, "bricks22.occ",
                                                         FORMAT, junk, CubitString(__FILE__) );
  ASSERT(s);

  //check that the two single volume bodys' attributes are exported as SINGLELUMP%
  Myfile.open ("bricks22.occ");
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

  assert (found == 2);

  // delete geometry
  GeometryQueryTool::instance()->delete_geometry();

  // import it again
  DLIList<RefEntity*> import_list;
  s = GeometryQueryTool::instance()->import_solid_model( "bricks2.occ",
                                                         FORMAT, NULL, CUBIT_FALSE,
                                                         CUBIT_TRUE, CUBIT_TRUE,
                                                         CUBIT_TRUE, CUBIT_TRUE,
                                                         CUBIT_TRUE, &import_list
                                                         );
  ASSERT(s);

  // check imported entity has tool data actuated by attributes
  DLIList<RefEntity*> body_entity_list;
  s = GeometryQueryTool::instance()->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
  int n_body = body_entity_list.size();
  body_entity_list.reset();
  assert (n_body == 2);

  export_list.clean_out();
  export_list.append(CAST_TO(body_entity_list.get_and_step(), Body));
  export_list.append(CAST_TO(body_entity_list.get_and_step(), Body));
  s = GeometryQueryTool::instance()->export_solid_model( export_list, "bricks23.occ",
                                                         FORMAT, junk, CubitString(__FILE__) );
  ASSERT(s);

  ASSERT(is_files_same("bricks2.occ", "bricks23.occ"));
  return 0;
}

CubitBoolean is_files_same(std::string filename1,
                           std::string filename2)
{
  FILE *fp1, *fp2;
  char ch1, ch2;
  CubitBoolean same = false;
  /* open first file */
  const char* file1 = filename1.c_str();
  const char* file2 = filename2.c_str();
  if((fp1 = fopen(file1, "rb"))==NULL) {
    printf("Cannot open first file.\n");
    exit(false);
  }

  /* open second file */
  if((fp2 = fopen(file2, "rb"))==NULL) {
    printf("Cannot open second file.\n");
    exit(false);
  } 

  /* compare the files */
  while(!feof(fp1)) {
    ch1 = fgetc(fp1);
    if(ferror(fp1)) {
      printf("Error reading first file.\n");
      exit(false);
    }
    ch2 = fgetc(fp2);
    if(ferror(fp2)) {
      printf("Error reading second file.\n");
      exit(false);
    }
    if(ch1 != ch2) {
      return same;
    }
  }

  if(fclose( fp1 ) == EOF) {
    printf("Error closing first file.\n");
    exit(false);
  }

  if(fclose( fp2 ) == EOF) {
    printf("Error closing second file.\n");
    exit(false);
  } 
  return true;
}

