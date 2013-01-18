#ifndef COMPOSITE_ENGINE_HPP
#define COMPOSITE_ENGINE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitString.hpp"
#include <map>
#include "DLIList.hpp"
#include "IntermediateGeomEngine.hpp"
#include "VGDefines.h"
#include "TopologyBridge.hpp"

class TBPoint;
class Curve;
class Surface;
class Lump;

class BodySM;
class ShellSM;
class LoopSM;
class CoEdgeSM;

class CompositePoint;
class CompositeCurve;
class CompositeSurface;
class CompositeLump;
class CompositeBody;
class CompositeCoEdge;
class CompositeLoop;
class CompositeShell;
class Body;

class GeometryEntity;
class TopologyBridge;
class CubitSimpleAttrib;
class CubitVector;
class CubitTransformMatrix;

class CompositeEngine : public IntermediateGeomEngine
{
	public:
    void get_all_curves_and_points(DLIList<TopologyBridge*> &tb_list,
                                   DLIList<Curve*> &curves,
                                   DLIList<TBPoint*> &points);
    bool is_composite(TBOwner *bridge_owner);
    bool is_composite(TopologyBridge *bridge);
    bool is_partition(TBOwner *bridge_owner);

    virtual void remove_imprint_attributes_after_modify
                                ( DLIList<BodySM*> &old_sms,
                                DLIList<BodySM*> &new_sms );
    virtual void push_named_attributes_to_curves_and_points
                     ( DLIList<TopologyBridge*> &tb_list, const char *name_in);
    virtual void push_imprint_attributes_before_modify
                     ( DLIList<BodySM*> &body_sms );
    virtual void attribute_after_imprinting(DLIList<TopologyBridge*> &tb_list,
                                            DLIList<Body*> &old_bodies);

    // This is a copy of the function in MergeTool with the difference that it
    // accepts a layer flag to dictate at which level the topology is traversed
    // (solid modeler level or virtual level).
    CubitBoolean about_spatially_equal( Curve *curve_1, Curve *curve_2,
                                                  CubitSense &relative_sense, 
                                                  double tolerance_factor,
                                                  int layer = TopologyBridge::MAX_TB_LAYER);
    virtual void remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges);

    virtual ~CompositeEngine();
	
    static CompositeEngine& instance();
    static void delete_instance();

    int level() const { return COMPOSITE_LAYER; }
		
    CubitStatus import_geometry( DLIList<TopologyBridge*>& imported_geometry );
      //- Recreate saved composite geometry and update the passed
      //- list.

    int is_hidden(TopologyBridge *tb);
    
    CubitStatus export_geometry( DLIList<TopologyBridge*>& geometry_to_save );
      //- Save composite geometry and update the passed list such
      //- that it contains the underlying real geometry.
		
    CompositeCurve* composite( Curve* survivor, 
                               Curve* dead,
                               TBPoint* keep_point = 0,
                               bool remove_partitions = false );
      //- Combine two curves.  The curves share both end vertices,
      //- the retained vertex can be specified in keep_point.
      
    CompositePoint* stitch_points( TBPoint* point1, TBPoint* point2 );
    CompositeCurve* stitch_curves( Curve* curve1, Curve* curve2 );
    CompositeSurface* stitch_surfaces( Surface* surf, Surface* surf2 );
    
    CompositeCurve* remove_point( TBPoint* dead_point,
                                  bool remove_partitions = false,
                                  Curve* survivor = 0 );
      //- Create a composite curve by removing the specified TBPoint.
                                    
                                    
    CompositeSurface* remove_curve( Curve* dead_curve,
                                    bool remove_partitions = false,
                                    Surface* survivor = 0 );
      //- Create a composite surface by removing the specified Curve
    
    CompositeLump* remove_surface( Surface* dead_surf,
                                   Surface* stitch_partner = 0,
                                   bool remove_partitions = false );
      //- Create a composite lump by removing the specified Surface.
      //- Note:  If stitch_partner is non-null, you need to call
      //-        combine_bodies on the owning bodies first.
                                   
    CompositeBody* combine_bodies( BodySM* body1, BodySM* body2 );
      //- Combine two BodySMs into one Composite Body.
    
    CubitStatus split_body( CompositeBody* body_to_split,
                            DLIList<BodySM*>& resulting_bodies );
      //- Split a CompositeBody into several Bodies 
      //- as permissible by common child topology and 
      //- existing real Bodies.
    
    CubitStatus restore_point( TBPoint* point );
      //- Split a composite curve by restoring the
      //- passed, previously removed TBPoint.
    
    CubitStatus restore_curve( Curve* curve );
      //- Split a composite surface by restoring the
      //- passed, previously removed Curve.
    
    CubitStatus restore_surface( Surface* surface, Surface*& stitch_parnter );
      //- Split a composite lump by restoring the
      //- passed, previously removed surface.
      
/*      
    static void fix_up_query_results( DLIList<TopologyBridge*>& list,
                                      bool keep_hidden_entities = false );
      //- Update TB-level query results obtained from calling 
      //- get_parents/children_virt such that the results are at the
      //- "composite level" of the TB graph.
*/    
    static CubitSimpleAttrib
    find_attribute_by_name( TopologyBridge* bridge, const CubitString name );     
      //- Find the attribute with the specified name.
      //- Caller is responsible for freeing the returned, heap-allocated
      //- CubitSimpleAttrib object.
      
    CompositePoint* replace_point( TBPoint* dead_point );
    CompositeCurve* replace_curve( Curve* dead_curve );
    CompositeSurface* replace_surface( Surface* surface );
    CompositeLump* replace_lump( Lump* lump );
    CompositeBody* replace_body( BodySM* body );
      //- These methods are public for use by other clases internal
      //- to the VG code.  
      
      /**************** To be called by PartitionEngine *****************/
      
    TBPoint* insert_point(CompositeCurve* curve, double u);
      //- Partition a composite curve.
      
    TBPoint* insert_point_curve( CompositeSurface* surface,
                               const CubitVector& position );
    
    CubitStatus insert_curve( DLIList<Surface*>& surfaces,
                              DLIList<CubitVector*>& polyline,
                              DLIList<Surface*>& new_surfaces,
                              DLIList<Curve*>& new_curves );

    void notify_deactivated (CompositeBody* body);
    void notify_deactivated (CompositeLump* volume);
    void notify_deactivated (CompositeSurface* surface);
    void notify_deactivated (CompositeCurve* curve);
    void notify_deactivated (CompositePoint* point);
    void clean_out_deactivated_geometry();

    CubitStatus restore_point_in_curve( TBPoint* point );
      //- Restore a point hidden by a curve, splitting the
      //- curve.
      
    CompositeCurve* restore_point_in_surface( TBPoint* point );
      //- Restore a point hidden by a surface, creating a 
      //- point-curve.


    CubitStatus translate( CompositeBody* body, const CubitVector& delta );
    CubitStatus rotate( CompositeBody* body, const CubitVector& axis, double degrees );
    CubitStatus scale( CompositeBody* body, const CubitVector& factors );
    CubitStatus reflect( CompositeBody* body, const CubitVector& axis );
    CubitStatus restore_transform( CompositeBody* body );
    
    CubitStatus translate( CompositeSurface* surface, const CubitVector& delta );
    CubitStatus rotate( CompositeSurface* surface, const CubitVector& axis, double degrees );
    CubitStatus scale( CompositeSurface* surface, const CubitVector& factors );
    CubitStatus reflect( CompositeSurface* surface, const CubitVector& axis );
    
    CubitStatus translate( CompositeCurve* curve, const CubitVector& delta );
    CubitStatus rotate( CompositeCurve* curve, const CubitVector& axis, double degrees );
    CubitStatus scale( CompositeCurve* curve, const CubitVector& factors );
    CubitStatus reflect( CompositeCurve* curve, const CubitVector& axis );

    CubitStatus notify_transform( TopologyBridge* bridge,
                                  const CubitTransformMatrix& xform );

    static void strip_attributes( TopologyBridge* bridge );
      // remove all attributes related to composite geometry
      // from passed bridge.

    void remove_attributes( DLIList<TopologyBridge*> &bridge_list );
      //remove Composite attributes off of topology bridges
    virtual void remove_modified(DLIList<Surface*> &all_surfs,
      DLIList<Curve*> &all_curves, DLIList<TBPoint*> &all_pts);

    void get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                               DLIList<TopologyBridge*> &tbs );
    
    
  protected:
	
    TBPoint* remove_composite( CompositePoint* point );
    Curve* remove_composite( CompositeCurve* composite );
    Surface* remove_composite( CompositeSurface* composite );
    Lump* remove_composite( CompositeLump* composite );
    BodySM* remove_composite( CompositeBody* composite );
    
    CubitStatus remove_partition_point( CompositePoint* point_owner );
    CubitStatus remove_partition_curves( CompositeCurve* curve_owner );
    
    CompositePoint* destroy_point_curve( CompositeCurve* curve );
    
    CompositeCurve* combine( CompositeCurve* curve1,
                             CompositeCurve* curve2,
                             CompositePoint* keep_point = 0,
                             bool remove_partitions = false );
    CubitStatus split( CompositeCurve* curve, int after_index,
                       Curve*& result1, Curve*& resutl2 );
                       
      // Shell modification helpers for composite lump creation
    CompositeShell* split_shell( CompositeShell* shell );
    void insert_nonmanifold_surfaces( DLIList<CompositeSurface*>& surfs,
                                      CompositeShell* shell1,
                                      CompositeShell* shell2 );
    CubitStatus inside_shell( CompositeShell* const shell,
                              CompositeSurface* const surf,
                              bool& result );
    
    CubitStatus create_composites( DLIList<BodySM*>& );
    CubitStatus create_composites( DLIList<Surface*>& );
    CubitStatus create_composites( DLIList<Curve*>& );
    CubitStatus create_composites( DLIList<TBPoint*>& );
      //- Methods called by import(..) to re-create
      //- composites from saved attributes.
      
    CubitStatus save( CompositePoint* );
    CubitStatus save( CompositeCurve* );
    CubitStatus save( CompositeSurface* );
    CubitStatus save( CompositeLump* );
    CubitStatus save( CompositeBody* );
    
    static void append_attrib( TopologyBridge* bridge, 
                               const CubitSimpleAttrib& attrib );
    
    
    static CubitStatus find_coedges( CompositeSurface* surface,
                                     CompositeCurve* curve,
                                     CompositePoint* point,
                                     CompositeCoEdge*& previous,
                                     CompositeCoEdge*& next );
    static CompositeCoEdge* find_next_point_coedge( 
                                     CompositeSurface* const surface,
                                     CoEdgeSM* const coedge,
                                     CompositePoint* point,
                                     DLIList<CompositeCoEdge*>& point_coedges );
    
    static CompositeEngine* instance_;
		
  private:
    
    CompositeEngine();
  
    CompositeSurface* split_surface( CompositeSurface* surf_to_split,
                                     CompositeLoop* loop_on_surf,
                                     CompositeLoop* new_loop );
    
    CompositeLump* split_lump( CompositeShell* shell_to_split );
      // If a shell had a previously-hidden surface inserted into
      // it due to the un-hiding of that surface, test if the
      // shell is split by the surface, and if so, split the shell
      // and corresponding lump.  
    void process_curves_after_imprint(Curve *att_bridge, 
                                                   Curve *other_bridge,
                                                   DLIList<BodySM*> &new_sms);
    void process_points_after_imprint(TBPoint *att_bridge, 
                                                   TBPoint *other_bridge,
                                                   DLIList<BodySM*> &new_sms);
    
    DLIList<TopologyBridge*> deactivatedList;
};

#endif
