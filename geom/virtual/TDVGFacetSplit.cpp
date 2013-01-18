#include "TDVGFacetSplit.hpp"
#include "ToolDataUser.hpp"

void TDVGFacetSplit::remove( ToolDataUser* user )
{
  delete user->remove_TD( &is_vg_facet_split );
}

