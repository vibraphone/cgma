#include "NoAcisModifyEngine.hpp"
#include "NoAcisQueryEngine.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitAttribManager.hpp"
#include "CGMApp.hpp"
#include "Body.hpp"
#include "CubitCompat.hpp"
#include <iostream>

extern int snippet();

int main (int argc, char **argv) 
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <acis_filename>" << std::endl;
    return 0;
  }
  
  CGMApp::instance()->startup(argc, argv);

    // declare acis engines
  AcisQueryEngine::instance();
  AcisModifyEngine::instance();
  

    // do something
  CGMApp::instance()->attrib_manager()->auto_flag(true);
  
  CubitStatus result = CubitCompat_import_solid_model(argv[1],"ACIS_SAT");
  if (CUBIT_SUCCESS != result) {
    std::cout << "Trouble opening file." << std::endl;
    return 1;
  }

  int dum = snippet();
  
  return dum;
}
