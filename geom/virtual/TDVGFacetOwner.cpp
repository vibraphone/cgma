#include "TDVGFacetOwner.hpp"
#include "ToolDataUser.hpp"

PartitionEntity* TDVGFacetOwner::get( ToolDataUser* user )
{
  ToolData* td = user->get_TD( &is_vg_facet_owner );
  return td ? dynamic_cast<TDVGFacetOwner*>(td)->owner() : 0;
}

void TDVGFacetOwner::set( ToolDataUser* user, PartitionEntity* owner )
{
  ToolData* td = user->get_TD( &is_vg_facet_owner );
  if( td )
  {
    dynamic_cast<TDVGFacetOwner*>(td)->mOwner = owner;
  }
  else
  {
    td = new TDVGFacetOwner( owner );
    user->add_TD( td );
  }
}

void TDVGFacetOwner::remove( ToolDataUser* user )
{
  delete user->remove_TD( &is_vg_facet_owner );
}

  
