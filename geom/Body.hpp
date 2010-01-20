/*-------------------------------------------------------------------------
 * Filename      : Body.hpp
 *
 * Purpose       : This class represents the topological entity, Body,
 *                 which is the highest level abstraction for a complete
 *                 geometric model.
 *
 * Special Notes : Body is a GroupingEntity in the TopologyEntity hierarchy.
 *                 It is also a RefEntity and get a lot of its functionality
 *                 and data from the RefEntity base class for meshing
 *                 purposes.
 *
 *                 The user is provided a handle to Bodies.
 *
 * Creator       : 
 *
 * Creation Date : 
 *
 * Owner         : 
 *-------------------------------------------------------------------------
 */

#ifndef BODY_HPP
#define BODY_HPP

#include "GroupingEntity.hpp"
#include "RefEntity.hpp"
#include "CoVolume.hpp"
#include "RefVolume.hpp"

class CubitTransformMatrix;

//! Body class.
class CUBIT_GEOM_EXPORT Body : public GroupingEntity,
             public RefEntity
{
public :

  /*- the factory is allowed to call the (private) constructors */
  friend class RefEntityFactory;

    /* constructors/destructors */

  virtual ~Body() ;
  /*- The destructor. */

    /* topology */
  
  //! Gets the dag type.
  DagType dag_type() const { return DagType::body_type(); }

  //! Gets the type info.
  const type_info& entity_type_info() const { return typeid(Body); }

  //! Gets the class name: "Body". 
  static const char* get_class_name()
     {
       return "Body";
     }

  //! Gets the class name: "Body". 
  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
  //! Gets the underlying BodySM pointer.
  BodySM* get_body_sm_ptr() const;

  //! Returns the bounding box of this body
  virtual CubitBox bounding_box();

  //! Return a CubitVector set to the centroid of this body
  virtual CubitVector center_point();

  //! Get certain mass props, cofg is center of gravity
  CubitBoolean get_mass_props(CubitVector& cofg);
       
  //! Determines whether a point is inside, outside, or on boundary
  //! of a volume.
  CubitPointContainment point_containment( CubitVector &point );
  
  //! Returns the volume of this body.
  virtual double measure();

  //! Query to see if this is a sheet body.
  CubitBoolean is_sheet_body();

#ifdef BOYD14
  CubitStatus get_transforms(CubitTransformMatrix &tfm);
  /*- returns the body transformation matrix
   */
#endif
  
#ifdef BOYD14
  CubitStatus transform_position(CubitVector const& position_vector,
                                 CubitVector & transformed_vector);
  /*R CubitStatus
   *R- CUBIT_SUCCESS/FAILURE
   *I position_vector
   *I- The input vector that needs to be transformed.
   *O transformed_vector
   *O- The transformed CubitVector -- sent in by reference.
   *- This function will transform the input position_vector using
   *- the transformation matrix associated with the underlying geometry
   *- of this Body.  If there is none, then CUBIT_FAILURE is returned
   *- and transformed_vector is left untouched.
   *- Note: The input vector *must* have 3 components.
   */
#endif

  //! Do a measure and api entity check.
  virtual int validate();
  
  //! Functions related to "measuring" the Body
  virtual CubitString measure_label();

  //! Sets the color.
  virtual void color(int value);
  //! Gets the color.
  virtual int color() const;

protected: 

  Body() ;
  /*- The default constructor.
   */

  Body(BodySM* OSMEPtr) ;
  /*- The constructor with a pointer to an other solid model entity.
   */

private:
   
#ifdef BOYD14
  void notify_model(CubitEventType eventType);
  /*R void
   *I eventType
   *I- The type of event that occurred 
   *- This function notifies the owning Model of the event that
   *- occured to this Body (see the EventType enum for the list of
   *- recognized events).
   */
#endif

  Body( const Body& );
  void operator=( const Body& );

} ;

#endif

