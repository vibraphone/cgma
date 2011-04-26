

#include "TestUtilities.hpp"
#include "TestConfig.h"

std::string data_file(char* filename)
{
  return DataDir + std::string("/") + std::string(filename);
}


bool cubit_box_identical(const CubitBox& box1, const CubitBox& box2, double tol,
    bool print_data)
{
  if(print_data)
  {
    printf("box 1 {(%f %f %f), (%f %f %f)}\n",
            box1.minimum().x(), 
            box1.minimum().y(), 
            box1.minimum().z(), 
            box1.maximum().x(), 
            box1.maximum().y(), 
            box1.maximum().z()
            );
    printf("box 2 {(%f %f %f), (%f %f %f)}\n",
            box2.minimum().x(), 
            box2.minimum().y(), 
            box2.minimum().z(), 
            box2.maximum().x(), 
            box2.maximum().y(), 
            box2.maximum().z()
            );
  }
  
  return box1.maximum().within_tolerance(box2.maximum(), tol) && 
         box1.minimum().within_tolerance(box2.minimum(), tol);
}


