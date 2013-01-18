//- Class:          CACompositeVG
//- Owner:          Tim Tautges
//- Description:    Cubit attribute for composite virtual geometry
//- Checked by:
//- Version:

#ifndef CA_COMPOSITE_VG_HPP
#define CA_COMPOSITE_VG_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

class Point;
class RefEntity;
class CubitSimpleAttrib;
class RefEntity;
class CubitVector;

class CACompositeVG: public CubitAttrib
{

private:

  int compositeId;
    //- the unique id of this composite entity

  DLIList<int> subEntityIds;
    //- the ids of the sub entities making up this composite

public:
  CACompositeVG(RefEntity*, const CubitSimpleAttrib&);
    //- construct a CACVG from a simple attribute

  virtual ~CACompositeVG() {};

  CubitStatus actuate();
  
  CubitStatus update();

  CubitStatus reset();

  CubitSimpleAttrib cubit_simple_attrib();
  
  CubitSimpleAttrib cubit_simple_attrib(CubitString);
  
  int int_attrib_type() {return CA_COMPOSITE_VG;}
    //- returns the enumerated attribute type

  int composite_id() {return compositeId;};
  void composite_id(int id) {compositeId = id;};
    //- get/set the composite id
  
  void check_child_cacvgs(RefEntity *new_entity);
    //- for the new entity passed in, check for and actuate CACVG's on
    //- any child entities
};

CubitAttrib* CACompositeVG_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif
