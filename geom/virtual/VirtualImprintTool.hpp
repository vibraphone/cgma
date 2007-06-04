//---------------------------------------------------------------------------
///
/// Class:          VirtualImprintTool
/// Description:    Performs the necessary tasks to imprint two surfaces
///                 and create new surfaces from the result.
/// Owner:          David R. White
/// Created:        7/11/2001
///
//---------------------------------------------------------------------------
#ifndef VIRTUALIMPRINTOOL_HPP
#define VIRTUALIMPRINTOOL_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
class RefFace;
class RefEdge;
class RefVertex;
class RefVolume;
class CubitVector;


template <class X> class DLIList;

class VirtualImprintTool 
{
protected:  
  VirtualImprintTool();
    ///
    /// Constructor.  This shouldn't
    /// get called since all of the functions are static.
    ///
  
  ~VirtualImprintTool();
    ///
    /// Destructor shouldn't do anything because
    /// it never allocates its own members
    ///
public:
  

  static CubitStatus virtual_imprint(RefFace *ref_face1,
                                     RefFace *ref_face2,
                                     DLIList <RefFace*> &results,
                                     double feature_size,
                                     CubitBoolean &curves_modified);
    ///
    /// Uses virtual geometry to imprint the two surfaces.
    ///

  static CubitStatus virtual_imprint(DLIList <RefFace*> &input_faces,
                                     DLIList <RefFace*> &surf_results,
                                     double feature_size,
                                     CubitBoolean &curves_modified);
    ///
    /// imprint multiple surfaces.
    ///

  
  static CubitStatus virtual_imprint(RefVolume *ref_volume1,
                                     RefVolume *ref_volume2,
                                     DLIList <RefVolume*> &vol_results,
                                     double feature_size,
                                     CubitBoolean &others_modified);
    ///
    /// Uses virtual geometry to imprint the two volumes.
    ///

  static CubitStatus virtual_imprint(DLIList <RefVolume*> &imprint_volumes,
                                     DLIList <RefVolume*> &vol_results,
                                     double feature_size,
                                     CubitBoolean &others_modified);
    ///
    /// Imprints the volumes in the list against each other.
    ///


private:
  static CubitBoolean useRealIntersection;
    ///
    /// Static flag storing preference on type of intersection for imprinting.
    ///

};
#endif

