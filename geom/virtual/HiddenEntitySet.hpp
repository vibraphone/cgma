//-------------------------------------------------------------------------
// Filename      : HiddenEntitySet.hpp
//
// Purpose       : A class to hold a list of hidden entites and act as
//                 their owner.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/14/02
//-------------------------------------------------------------------------

#ifndef HIDDEN_ENTITY_SET
#define HIDDEN_ENTITY_SET

#include "TBOwner.hpp"
#include <vector>

class BodySM;
class Lump;
class Surface;
class Curve;
class TBPoint;

class HiddenEntitySet : public TBOwner
{

  public:
  
    HiddenEntitySet( TBOwner* owner ) : myOwner(owner) {}
  
    virtual ~HiddenEntitySet();
  
    CubitStatus hide( TopologyBridge* bridge );
    CubitStatus restore( TopologyBridge* bridge );

    void hidden_surfaces( DLIList<Surface*>& );
    void hidden_coedges(DLIList<CoEdgeSM*>& );
    void hidden_curves( DLIList<Curve*>& );
    void hidden_points( DLIList<TBPoint*>& );
    
    CubitStatus remove_bridge( TopologyBridge* bridge );
    CubitStatus swap_bridge( TopologyBridge* old_tb,
                             TopologyBridge* new_tb,
                             bool reversed );
    void notify_reversed( TopologyBridge* );
    
    CubitStatus merge( HiddenEntitySet* other );
    
    TBOwner* owner() const
      { return myOwner; }
      
    void print_debug_info( const char* prefix ) const;
    
    bool is_empty() const
      { return hiddenList.empty(); }
    
  private:
  
    TBOwner* myOwner;
    std::vector<TopologyBridge*> hiddenList;
};

#endif
