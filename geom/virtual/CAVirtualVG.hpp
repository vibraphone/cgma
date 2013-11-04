//- Class:          CAVirtualVG
//- Owner:          Tim Tautges
//- Description:    Cubit attribute for virtual geometry entities
//-                 (holds both virtual vertices and virtual curves)
//- Checked by:
//- Version:
//-
//- This attribute stores:
//- numVV: number of virtual vertices
//- numVC: number of virtual curves
//- vgUIDs[0..numVV-1]: UID's of virtual vertices
//- vgUIDs[numVV..(end)]: for each curve, UID's of start, end vertices, then
//-                 UID of virtual curve
//- posVector[0..numVV-1]: position of virtual vertices
//- posVector[numVV..(end)]: positions of vpoints defining vcurve (not including ends)
//- numVCPoints: number of virtual points for each virtual curve

#ifndef CA_VIRTUAL_VG_HPP
#define CA_VIRTUAL_VG_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

class Point;
class RefEntity;
class CubitSimpleAttrib;
//class ParasitePoint;
//class ParasiteCurve;

class CAVirtualVG: public CubitAttrib
{

friend class CAPartitionVG;

protected:
  int numVV, numVC;
    //- the number of virtual points and curves contained in this attr

  DLIList<int> vgUIDs;
    //- unique ids of virtual points and curves contained in this attr

  DLIList<CubitVector*> posVector;
    //- position vectors for virtual curves

  DLIList<int> numVCPoints;
    //- for each virtual curve, the number of virtual points in that curve

//  void add_vvertex(ParasitePoint *vpoint);
    //- adds data for the RefVertex defined by this vpoint to this CA

//  void add_vcurve(ParasiteCurve *vcurve);
    //- adds data for this vcurve to this CA

public:

  CAVirtualVG(RefEntity*, const CubitSimpleAttrib&);
    //- construct a CAPVG from a simple attribute

  virtual ~CAVirtualVG() {};

  CubitStatus actuate();
  
  CubitStatus update();

  CubitStatus reset();
    //- reset info; called from CAU and also from update!

  CubitSimpleAttrib cubit_simple_attrib();
  
  CubitSimpleAttrib cubit_simple_attrib(CubitString);
  
  int int_attrib_type() {return CA_VIRTUAL_VG;}
    //- returns the enumerated attribute type

};

CubitAttrib* CAVirtualVG_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif
