#ifndef OTHER_REF_ENTITY_HPP
#define OTHER_REF_ENTITY_HPP

#include "RefEntity.hpp"

class OtherEntity;


class OtherRefEntity : public RefEntity
{
public:
  OtherRefEntity();
  virtual ~OtherRefEntity();
  OtherEntity* get_other_entity_ptr(){ return otherEntityPtr; }
  OtherEntity const* get_other_entity_ptr() const { return otherEntityPtr; }
    //R OtherEntity*
    //R- A pointer to a other entity  associated with the 
    //R- grouping entity.
    //- This function returns a pointer to the other solid model entity 
    //- associated with this object.
  
  CubitStatus set_other_entity_ptr(OtherEntity* other_entity_ptr);
  
protected:
  OtherEntity* otherEntityPtr;
protected:
  
};

#endif

