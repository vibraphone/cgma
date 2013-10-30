//- Class:          CAEntitySense
//- Owner:          Corey Ernst
//- Description:    Saves out bridge sense of RefEntity
//- Checked by:
//- Version:

#ifndef CA_ENTITY_SENSE_HPP
#define CA_ENTIYT_SENSE_HPP

#include <typeinfo>
#if !defined(WIN32)
using std::type_info;
#endif

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntitySense: public CubitAttrib
{

private:

  CubitSense entitySense;
  
public:
  CAEntitySense(RefEntity* = NULL);

  CAEntitySense(RefEntity*, const CubitSimpleAttrib&);
    //- make a CAG from a simple attribute

  void initialize();
    //- initialize random number generator for this attribute

  virtual ~CAEntitySense();

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAEntitySense);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();
    //- actuate this attribute
  
  CubitStatus update();
    //- update this attribute

  CubitStatus reset() {return CUBIT_SUCCESS;}
    //- reset this attribute

  CubitSimpleAttrib cubit_simple_attrib();
    //- return a simple attribute with this CA's data  
  
  int int_attrib_type() {return CA_ENTITY_SENSE;}
    //- returns the enumerated attribute type

  void print();

};

CubitAttrib* CAEntitySense_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

