#ifndef TD_VG_FACET_SPLIT_HPP
#define TD_VG_FACET_SPLIT_HPP

#include "ToolData.hpp"

class CubitFacetData;
class ToolDataUser;

class TDVGFacetSplit : public ToolData
{
  public:
  
    TDVGFacetSplit( CubitFacetData* split_from )
      : splitFrom( split_from ) {}
    
    static int is_vg_facet_split( const ToolData* td )
      { return !!dynamic_cast<const TDVGFacetSplit*>(td); }
    
    CubitFacetData* split_from( ) const
      { return splitFrom; }

    static void remove( ToolDataUser* user );
      
  private:
    
    CubitFacetData* const splitFrom;
};

#endif
