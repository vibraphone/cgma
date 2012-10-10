
#include "CubitCoordinateSystem.hpp"

#include <cassert>

#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "CubitObserver.hpp"

//! default constructor, initializes to rectangular type with an identity transform
CubitCoordinateSystem::CubitCoordinateSystem(int id)
 : mType(Rectangular), mMatrix(CubitTransformMatrix()), mReference(0)
{
    this->set_id(id);
}

CubitCoordinateSystem::~CubitCoordinateSystem()
{
  if(mReference)
    mReference->mUses.remove(this);
  
  // if we have uses, make them reference this reference
  while(mUses.size())
  {
    CubitCoordinateSystem* sys = mUses.get();
    sys->set_transform(sys->mMatrix * this->mMatrix);
    sys->mReference = this->mReference;
    if(mReference)
      mReference->mUses.append(sys);
    mUses.remove();
  }

  this->remove_from_observers();

  CubitObserver::notify_static_observers(this, COORDINATE_SYSTEM_DELETED);
}

CubitBox CubitCoordinateSystem::bounding_box()
{
    assert(0); // this function is only here b/c CubitEntity requires it
    
    CubitBox temp;
    return temp;
}

const type_info& CubitCoordinateSystem::entity_type_info() const
{
    return typeid(CubitCoordinateSystem);
}

const char* CubitCoordinateSystem::class_name() const
{
    return "Coordinate System";
}

CubitCoordinateSystem::Type CubitCoordinateSystem::get_type() const
{
  return mType;
}

void CubitCoordinateSystem::set_type(CubitCoordinateSystem::Type type)
{
  mType = type;
  CubitObserver::notify_static_observers(this, COORDINATE_SYSTEM_MODIFIED);
}


//! return the transform of this coordinate system
const CubitTransformMatrix& CubitCoordinateSystem::get_transform() const
{
  return mMatrix;
}

void CubitCoordinateSystem::set_transform(const CubitTransformMatrix& mat)
{
  mMatrix = mat;
  CubitObserver::notify_static_observers(this, COORDINATE_SYSTEM_MODIFIED);
}


//! return the transform of this coordinate system 
//! with respect to the global coordinate system
CubitTransformMatrix CubitCoordinateSystem::get_concatenated_transform() const
{
  CubitTransformMatrix mat;
  mat = mMatrix;
  CubitCoordinateSystem* ref;
  for(ref = mReference; ref != 0; ref = ref->mReference)
    mat = mat * ref->mMatrix;
  return mat;
}


