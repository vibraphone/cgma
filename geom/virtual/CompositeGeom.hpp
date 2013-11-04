//-------------------------------------------------------------------------
// Filename      : CompositeGeom.hpp
//
// Purpose       : Object used by CompositeSurface and CompositeCurve
//                 to manage underlying entities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------

#ifndef COMPOSITE_GEOM_HPP
#define COMPOSITE_GEOM_HPP

#include "DLIList.hpp"
#include "VGArray.hpp"
#include "CubitBox.hpp"

class GeometryEntity;
class CubitSimpleAttrib;
class CompositeAttrib;
class TopologyBridge;

struct CompositeEntry
{
  public:
  
    GeometryEntity* entity; // underlying entity
    CubitSense sense;       // sense relative to composite
    double measure;         // cached value of entity->measure()
    CubitBox bbox;          // cached value of entity->bounding_box()
    double dist_sqr;        // temporary data used in closest_point() calc.
    
    inline CompositeEntry()
      : entity(0), sense(CUBIT_UNKNOWN), measure(0.0), dist_sqr(0.0) {}
};


class CompositeGeom
{

  public:
  
    CompositeGeom( int size = 0 );
    
    ~CompositeGeom();
    
    int num_entities() const
      { return entityList.size(); }
    
    int index_of( TopologyBridge* geom_ptr ) const;
    inline GeometryEntity* entity( int index ) const;
    
    inline CubitSense sense( int index ) const;
    inline double measure( int index );
    void reverse_sense( int index );
    
    CubitStatus insert( int index, GeometryEntity* geom_ptr, CubitSense sense );
    
    CubitStatus append( GeometryEntity* geom_ptr, CubitSense sense )
      { return insert( entityList.size(), geom_ptr, sense ); }
    
    CubitStatus remove( int index, bool dead );
    CompositeGeom* split( int index );
    CompositeGeom* split( VGArray<int>& index_array );
    CubitStatus swap( int index, GeometryEntity* new_geom );
    CubitStatus reverse();
    CubitStatus reverse_order();
    CubitStatus reverse_rel_senses();
    
    CubitStatus merge( CompositeGeom& dead, bool prepend = false );
    
    CubitBox bounding_box();
    double measure();
    
    int next_box_within_dist( double dist_squared );
    int closest_box( const CubitVector& position );
    
    void update_cached_data()
      { needToUpdateBbox = true; needToUpdateMeasure = true; }
    
    void add_attribute( const CubitSimpleAttrib& csa );
    void rem_attribute( const CubitSimpleAttrib& csa );
    void rem_all_attributes();
    void get_attributes( DLIList<CubitSimpleAttrib>& list );
    void get_attributes( const char* name,
                         DLIList<CubitSimpleAttrib>& list );
    
    void print_debug_info( const char* line_prefix = 0 );
    
    void read_attributes( GeometryEntity* entity = 0 );
    void write_attributes( GeometryEntity* entity = 0 );
    
  private:
  
    static void clean_up_attribs( GeometryEntity* ent );

      // these have no implementation, just private delcarations
      // to prevent the compiler from generating default implementations
    CompositeGeom& operator=(const CompositeGeom&);
    CompositeGeom(const CompositeGeom&);
  
    void update_data_bbox();
    void update_data_measure();
    
    VGArray<CompositeEntry> entityList;
    int currentIndex;
    int firstIndex;
    
    bool needToUpdateBbox;
    bool needToUpdateMeasure;
    
    CompositeAttrib* listHead;
};

    
inline GeometryEntity* CompositeGeom::entity( int index ) const
{ 
  return entityList[index].entity;
}
    
inline CubitSense CompositeGeom::sense( int index ) const
{ 
  return entityList[index].sense;
}

inline double CompositeGeom::measure( int index ) 
{ 
  if( needToUpdateMeasure ) 
    update_data_measure();
  return entityList[index].measure;
}
    
#endif

     
    
    
    
