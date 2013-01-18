//-------------------------------------------------------------------------
// Class:       TDInterpNode
// Description: Support for the TetFacetorTool.  Maintains interpolation data

// Author:      sjowen
// Date:        8/13/2003
//-------------------------------------------------------------------------

#ifndef TD_INTERP_NODE_HPP
#define TD_INTERP_NODE_HPP

#include "ToolData.hpp"
#include "CastTo.hpp"


class TDInterpNode : public virtual ToolData
{
private:

  double mValue;
  double mWeight;

public:

  TDInterpNode( double val )
    { mValue = val; mWeight = 0.0; }
  virtual ~TDInterpNode(){};
   //-constructor and destructor

  static int is_interpnode(const ToolData* td)
     {return (CAST_TO(td, const TDInterpNode) != NULL);}

  void value( double val )
    { mValue = val; }
  double value()
    { return mValue; }

  void weight( double wgt )
    { mWeight = wgt; }
  double weight()
    { return mWeight; }
};



#endif // TD_INTERP_NODE_HPP

