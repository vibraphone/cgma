#include "cgm_test.hpp"
#include "CubitMessage.hpp"

int TEST(int argc, char** argv);

int main(int argc, char** argv)
{
  start_cgm(argc,argv);
  int result = TEST(argc,argv);
  int ret_val = ( CubitMessage::instance()->error_count() );
  
  static std::string cgm_port = argv[2];
  if(argc > 2 && cgm_port == "occ" )
    CubitMessage::instance()->reset_error_count(ret_val -5); 
  
  end_cgm(argc,argv);
  return result;
}
