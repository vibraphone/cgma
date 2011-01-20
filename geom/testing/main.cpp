#include "cgm_test.hpp"

int TEST(int argc, char** argv);

int main(int argc, char** argv)
{
  start_cgm(argc,argv);
  int result = TEST(argc,argv);
  end_cgm(argc,argv);
  return result;
}
