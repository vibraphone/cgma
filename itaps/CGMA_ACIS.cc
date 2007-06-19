#include "CGMA_ACIS.h"

class AcisModifyEngine 
{
public:
  AcisModifyEngine();
  static AcisModifyEngine *instance_;
};


class AcisQueryEngine 
{
public:
  AcisQueryEngine();
  static AcisQueryEngine *instance_;
};

class dummym
{
  char dummyvar[4096];
};

int CGMA_initialize_ACIS()
{
  have_initialized_acis = true;
  AcisQueryEngine ::instance_ = new (new dummym) AcisQueryEngine;
  AcisModifyEngine::instance_ = new (new dummym) AcisModifyEngine;
  return 0;
}

  
