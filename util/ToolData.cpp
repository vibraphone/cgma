// file ToolData.cc

#include "ToolData.hpp"

#include "CastTo.hpp"

ToolData::~ToolData() {}

// note virtual functions often can't be inlined, and the "is_xxx"
// functions aren't inlined in a lot of cases where their address is
// taken.

ToolData* ToolData::propogate(ToolDataUser* new_td_user)
{
  return NULL;
}

ToolData* ToolData::merge(ToolDataUser* other_td_user)
{
  return NULL;
}
