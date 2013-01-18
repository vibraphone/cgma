#ifndef PARTITION_SHELL_HPP
#define PARTITION_SHELL_HPP

#include "ShellSM.hpp"
#include "PartitionCoSurf.hpp"

class PartitionSurface;
class PartitionLump;
class CubitVector;

class PartitionShell : public ShellSM
{
friend class PartitionLump;

public:
  
  PartitionShell( );
  virtual ~PartitionShell();
  
  PartitionLump* get_lump() const;
  
  PartitionCoSurf* next_co_surface( const PartitionCoSurf* prev = 0 ) const;
  
  CubitStatus add( PartitionCoSurf* cosurf );
  CubitStatus remove( PartitionCoSurf* cosurf );
  
  PartitionCoSurf* add( PartitionSurface* surf, CubitSense sense );
    // create a CoSurf
  PartitionCoSurf* find_first( const PartitionSurface* surface ) const;
    // find first CoSurf with the passed surface
  PartitionCoSurf* find_next( const PartitionCoSurf* cosurf ) const;
    // find next CoSurf with the same surface
  CubitSense find_sense( const PartitionSurface* surface ) const;
    // returns CUBIT_UNKNOWN if multiple CoSurfs
  void remove_all_surfaces( DLIList<PartitionSurface*>* removed = 0 );
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  int layer() const;
  GeometryQueryEngine* get_geometry_query_engine() const;
  
  void append_simple_attribute_virt( const CubitSimpleAttrib&  );
  void remove_simple_attribute_virt( const CubitSimpleAttrib&  );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString&,
                                    DLIList<CubitSimpleAttrib>& );
  
  void print_debug_info( const char* prefix = 0 ) const;
  
  CubitPointContainment point_containment( const CubitVector& pt );
  
  bool is_nonmanifold( PartitionSurface* surface ) const;
  
  CubitStatus mass_properties( CubitVector& centroid, double& volume );
  
private:

  PartitionLump* myLump;
  PartitionShell* lumpNext;
  
  PartitionCoSurf* firstCoSurf;
};

inline PartitionLump* PartitionShell::get_lump() const
  { return myLump; }

inline PartitionCoSurf* 
PartitionShell::next_co_surface( const PartitionCoSurf* prev ) const
  { return !prev ? firstCoSurf : prev->myShell == this ? prev->shellNext : 0; }

#endif
