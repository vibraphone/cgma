//-------------------------------------------------------------------------
// Filename      : ChollaSkinTool.hpp 
//
// Purpose       : For a list of facets, finds the exterior edge facets.
//
// Special Notes : 
//
// Creator       : Steven Owen
//
// Date          : 4/27/01
//
// Owner         : Steven Owen
//-------------------------------------------------------------------------

#ifndef ChollaSkinTool_HPP
#define ChollaSkinTool_HPP

#include "DLIList.hpp"
class FacetEntity;
class ChollaSurface;
class ChollaCurve;
class ChollaPoint;

class ChollaSkinTool
{
public:
  ChollaSkinTool();
   //- Constructor

  ~ChollaSkinTool();
   //- Destructor

  CubitStatus skin_2d(DLIList<FacetEntity*> &facet_list,
                      ChollaSurface *&facet_surface_mesh_ptr );
    //- returns the bounding list of edges for the given set
    //- of mesh entities.  The mesh list can contain, quads and/or tris

  CubitStatus skin_1d(DLIList<FacetEntity*> &facet_list,
                      ChollaCurve *&facet_curve_mesh_ptr );
    //- returns the bounding list of nodes for the given set
    //- of mesh entities.  The mesh list can contain edges

private:

};


#endif

