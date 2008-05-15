//-------------------------------------------------------------------------
// Filename      : PartitionBody.hpp
//
// Purpose       : BodySM implementation for partition geometry
//
// Special Notes : Catches transforms
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/13/03
//-------------------------------------------------------------------------

#ifndef PARTITION_BODY_HPP
#define PARTITION_BODY_HPP

#include "BodySM.hpp"
#include "PartitionEntity.hpp"

class CubitFacetPoint;
class SegmentedCurve;
class PartitionPoint;

class PartitionBody : public BodySM, public PartitionEntity
{
  public:
  
    PartitionBody( BodySM* );
    
    ~PartitionBody();
    

    /************** Functions from BodySM **************/
/*    
    CubitStatus move( const CubitVector& offset );
    
    CubitStatus rotate( const CubitVector& axis, double angle );
    
    CubitStatus scale( double factor );
    
    CubitStatus scale( const CubitVector& factors );
    
    CubitStatus reflect( const CubitVector& axis );
    
    CubitStatus restore();
    
    CubitStatus reverse();
*/    
    CubitStatus get_transforms( CubitTransformMatrix& );
    
    CubitStatus mass_properties( CubitVector& centroid, double& volume );
    
    CubitPointContainment point_containment( const CubitVector& pos );
    
 
    /************** Functions from TopologyBridge **************/
    
    void append_simple_attribute_virt( CubitSimpleAttrib* );
    
    void remove_simple_attribute_virt( CubitSimpleAttrib* );
    
    void remove_all_simple_attribute_virt();
    
    CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib*>& );
    
    CubitStatus get_simple_attribute( const CubitString& name,   
                                      DLIList<CubitSimpleAttrib*>& );
    
    virtual int layer() const;
    
    void get_parents_virt( DLIList<TopologyBridge*>& );
    
    void get_children_virt( DLIList<TopologyBridge*>& );
    
    GeometryQueryEngine* get_geometry_query_engine() const;
    
    
     /************** Functions from PartitionEntity **************/
   
    void reverse_sense();
    
    void notify_split( FacetEntity*, FacetEntity* );
    
    CubitStatus save( CubitSimpleAttrib& );
    
    CubitBox bounding_box() const;
    
    
     /************** Local Functions **************/
    
    BodySM* real_body() const;
    
    void add( SubEntitySet& );
    
    void remove( SubEntitySet& );
    
    bool has_children() const { return !!childList; } 
    
    void destroy_all_children();
    
    void transform(const CubitTransformMatrix&) {}
    
    void get_all_children( DLIList<PartitionEntity*>& list );

  protected:
    
  private:
  
    SubEntitySet* childList;
};

#endif

    
