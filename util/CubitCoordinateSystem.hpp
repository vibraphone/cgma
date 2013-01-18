
#ifndef CUBIT_COORDINATE_SYSTEM_HPP
#define CUBIT_COORDINATE_SYSTEM_HPP


#include "CubitTransformMatrix.hpp"
#include "DLIList.hpp"
#include "CubitUtilConfigure.h"
#include "CubitEntity.hpp"
#include "CubitObservable.hpp"

class CUBIT_UTIL_EXPORT CubitCoordinateSystem : public CubitEntity,
                                                public CubitObservable
{
  public:
    
    //! coordinate systems are of the types:
    enum Type { Rectangular, Cylindrical, Spherical };

    //! default constructor, initializes to rectangular type with an identity transform
    CubitCoordinateSystem(int id);
    
    //! destructor, if this system is concatenated to others, their matrices will
    //! by multiplied by this transform
    ~CubitCoordinateSystem();

    //! get the id
    //int id() const { return mId; } use CubitEntity's ID
    virtual CubitBox bounding_box();
    
    virtual const type_info&  entity_type_info() const;
    virtual const char* class_name() const;
    //- return class name string.

    virtual CubitString entity_name() const;
    //- return entity name

    //! get and set the type of this coordinate system
    Type get_type() const;
    void set_type(Type);
    
    //! return the transform of this coordinate system
    const CubitTransformMatrix& get_transform() const;

    //! set the transform of this coordinate system
    void set_transform(const CubitTransformMatrix& mat);

    //! return the transform of this coordinate system 
    //! with reference to the global coordinate system
    CubitTransformMatrix get_concatenated_transform() const;
  
    
  protected:
    // the id of this system
    //int mId; now use CubitEntity's ID
    // the type of this system
    Type mType;
    // the transform that represents this coordinate system's transformation
    CubitTransformMatrix mMatrix;
    // this coordinate system may reference other coordinate systems
    CubitCoordinateSystem* mReference;
    // this coordiante system may be referenced by multiple other systems
    DLIList<CubitCoordinateSystem*> mUses;

};

#endif //CUBIT_COORDINATE_SYSTEM_HPP

