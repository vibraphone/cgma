#include "GlobalCommandFeedback.hpp"

namespace GlobalCommandFeedback
{
  ManagedPtrVector<CommandFeedback>* instance_;
  
  void create()
  { instance_ = new ManagedPtrVector<CommandFeedback>; }
    
  void destroy()
  {
    delete instance_;
    instance_ = NULL;
  }
  
  ManagedPtrVector<CommandFeedback>& instance()
  { return *instance_; }
}
