
#include "stdio.h"
#include "TestUtilities.hpp"

#include "GeometryQueryTool.hpp"

int ReadIgesFile(int argc, char* argv[])
{

  std::string file = data_file("brick.iges");

  CubitStatus stat = GeometryQueryTool::instance()->import_solid_model(file.c_str(), "IGES");

  if(stat != CUBIT_SUCCESS)
  {
    printf("failed to read iges file %s\n", file.c_str());
    return 1;
  }

  return 0;
}

