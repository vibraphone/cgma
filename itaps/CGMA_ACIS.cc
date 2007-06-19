#include "CGMA_ACIS.h"
#include <new>

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
  if (!AcisQueryEngine::instance_) {
    AcisQueryEngine::instance_ = reinterpret_cast<AcisQueryEngine*>(new dummym);
    new (AcisQueryEngine::instance_) AcisQueryEngine;
  }
  if (!AcisModifyEngine::instance_) {
    AcisModifyEngine::instance_ = reinterpret_cast<AcisModifyEngine*>(new dummym);
    new (AcisModifyEngine::instance_) AcisModifyEngine;
  }
  return 0;
}

  
