#ifndef TD_VG_FACET_OWNER_HPP
#define TD_VG_FACET_OWNER_HPP

#include "ToolData.hpp"

class PartitionEntity;
class ToolDataUser;

class TDVGFacetOwner : public ToolData
{

  public:
  
    TDVGFacetOwner( PartitionEntity* owner )
      : mOwner( owner ) {}
    
    static int is_vg_facet_owner( const ToolData* td )
      { return !! dynamic_cast<const TDVGFacetOwner*>(td); }
    
    PartitionEntity* owner() const
      { return mOwner; }
    
    static PartitionEntity* get( ToolDataUser* user );
    static void set( ToolDataUser* user, PartitionEntity* owner );
    static void remove( ToolDataUser* user );
  
  private:
  
    PartitionEntity* mOwner;
};

#endif
