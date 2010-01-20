

#include "cgm_test.hpp"

#include "AppUtil.hpp"
#include "TextProgressTool.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"

#ifdef CGM_SMLIB
#include "SMLibQueryEngine.hpp"
#include "SMLibModifyEngine.hpp"
#endif

#ifdef CGM_ACIS
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#endif


static std::string cgm_port;

void start_cgm(int argc, char** argv)
{
  for(int i=0; i<argc; i++)
  {
    // get library name
    if(std::string(argv[i]) == "-E" && i+1 < argc)
    {
      cgm_port = argv[i+1];
    }
  }


  // init util
  AppUtil::instance()->startup(argc, argv);
  AppUtil::instance()->progress_tool(new TextProgressTool());

  // init cgm
  CGMApp::instance()->startup(argc, argv);

  GeometryQueryEngine* query_engine = NULL;
  GeometryModifyEngine* modify_engine = NULL;
  printf("using %s engine\n", cgm_port.c_str());

  // find the engine
  if(0)
  {
    // nothing here
  }
#ifdef CGM_ACIS
  else if(cgm_port == "acis")
  {
    query_engine = AcisQueryEngine::instance();
    modify_engine = AcisModifyEngine::instance();
  }
#endif
#ifdef CGM_SMLIB
  else if(cgm_port == "smlib")
  {
    query_engine = SMLibQueryEngine::instance();
    modify_engine = SMLibModifyEngine::instance();
  }
#endif

  if(query_engine == NULL || modify_engine == NULL)
  {
    printf("no engine set with -E flag\n");
    exit(1);
  }
}

void end_cgm(int argc, char** argv)
{
  CGMApp::instance()->shutdown();
  CGMApp::delete_instance();

  // close down util module
  AppUtil::instance()->shutdown();
  AppUtil::delete_instance();
}

