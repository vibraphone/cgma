//-------------------------------------------------------------------------
// Filename      : PartitionEngine.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/23/02
//-------------------------------------------------------------------------

#ifndef PARTITION_ENGINE_HPP
#define PARTITION_ENGINE_HPP

#include "CubitDefines.h"
#include "VGDefines.h"
#include <map>
#include "IntermediateGeomEngine.hpp"
#include "DLIList.hpp"

class TopologyBridge;
class TBPoint;
class Curve;
class Surface;
class Lump;
class LoopSM;
class BodySM;

class SubSurface;
class SubCurve;
class SegmentedCurve;

class PartitionEntity;
class PartitionBody;
class PartitionLump;
class PartitionShell;
class PartitionCoSurf;
class PartitionSurface;
class PartitionCurve;
class PartitionPoint;
class PartitionLoop;
class PartitionCoEdge;
class PartPTCurve;

class PST_Point;
class SubEntitySet;
class FacetProjectTool;
class CubitSimpleAttrib;

class CubitPoint;
class CubitFacet;
class CubitPointData;
class CubitFacetData;
class CubitFacetEdgeData;
class Body;

class CubitVector;
class CubitTransformMatrix;
template <class X> class DLIList;

class PartitionEngine : public IntermediateGeomEngine
{

public:
  bool is_partition(TBOwner *bridge_owner);
  bool is_composite(TBOwner *bridge_owner);
  bool is_composite(TopologyBridge *bridge);
  virtual void remove_imprint_attributes_after_modify
                                ( DLIList<BodySM*> &old_sms,
                                DLIList<BodySM*> &new_sms ){};
  virtual void push_imprint_attributes_before_modify
                     ( DLIList<BodySM*> &body_sms ){};
  virtual void push_named_attributes_to_curves_and_points
                     ( DLIList<TopologyBridge*> &tb_list, const char *name_in ){};
  virtual void attribute_after_imprinting(DLIList<TopologyBridge*> &tb_list,
                                          DLIList<Body*> &old_bodies){};
  virtual void remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges){};



  virtual ~PartitionEngine();

  /** Get singleton instance */
  static PartitionEngine& instance();
  static void delete_instance();

  int level() const { return SUBCOMP_PARTITION_LAYER; }

  void notify_deactivated (PartitionBody* body);
  void notify_deactivated( PartitionLump* vol);
  void notify_deactivated (PartitionSurface* surface);
  void notify_deactivated (PartitionCurve* curve);
  void notify_deactivated (PartitionPoint* point);
  void add_to_deactivated_list(PartitionBody* body);
  void add_to_deactivated_list(PartitionLump* vol);
  void add_to_deactivated_list(PartitionSurface* sur);
  void add_to_deactivated_list(PartitionCurve* cur);
  void add_to_deactivated_list(PartitionPoint* pt);
  /** Restore partitions on passed geometry
   *
   * Reads partition geometry saved as attributes on the passed
   * geometry, restores the partitions and udpates the passed
   * list
   */
  CubitStatus import_geometry( DLIList<TopologyBridge*>& imported_geometry );

  /** Save partition geometry.
   *
   * Saves any partition geometry in the passed list (or on children
   * of entities int the passed list) as attributes on the underlying
   * geometry, and update the passed list to contain the underlying
   * geometry.
   */
  CubitStatus export_geometry( DLIList<TopologyBridge*>& geometry_to_save );

  /** Partition a curve
   *
   * @param curve The curve to partition
   * @param u     The location at which to split the curve.
   * @return      The new point created at the split location.
   */
  TBPoint* insert_point( Curve* curve, double u );
  
  /** Undo curve partitioning or point-curve creation.
   * 
   * @param point The point to remove.
   * @param dead_curves If non-null, polulated with pointers to the
   *                    destroyed PartitionCurves (pointers are by
   *                    definition stale pointers.)
   * @return The surviving curve, or null if a point-curve was removed..
   */
  Curve* remove_point( PartitionPoint* point, 
                       PartitionCurve* dead_curves[2] = 0 );
  
  /** Single-surface interface to insert_curve(DLIList<Surface*>&,...)
   *
   * @see insert_curve( DLIList<Surface*>&, DLIList<CubitVector*>&, DLIList<Surface*>,& DLIList<Curve*>& ) 
   */                   
  Surface* insert_curve( Surface* surface, 
                         DLIList<CubitVector*>& segment_points,
                         DLIList<Curve*>& new_curves,
                         const double *tolerance_length = NULL);
  
  /** Partition surfaces
   * 
   * Split surfaces by projecting the passed polyline onto the patch
   * of surfaces.
   * @param surfaces     The list of surfaces to split.
   * @param polyline     The curve to imprint the surfaces with.
   * @param new_surfaces Populated with new surfaces created during partitioning.
   * @param polyline_curves Popolated with new curves.
   */
  CubitStatus insert_curve( DLIList<Surface*>& surfaces,
                            DLIList<CubitVector*>& polyline,
                            DLIList<Surface*>& new_surfaces,
                            DLIList<Curve*>& polyline_curves,
                            const double *tolerance_length = NULL,
                            DLIList<Surface*>* surfs_to_reverse = NULL);
  
  /** Imprint a point on a surface, creating a zero-length (point) curve.
   *
   * Imprint a point on a surface, creating a zero-length (point) curve.
   * @param surface  The surface on which to create the point-curve
   * @param position The location on the surface at which to create the point.
   * @return         The point created (and owned by the point-curve).  The
   *                 point curve will be the only parent curve of this point.
   */
  TBPoint* insert_point_curve( Surface* surf,
                             const CubitVector& position,
                             Surface *&partitioned_surf );

  /** Partition a lump.
   *
   * Create a new surface in the lump as a copy of the passed surface.
   * This function may partition the lump if the surface copy divides
   * a shell of the lump into two distinct regions.
   * @param surface_to_copy  The geometry from which to create the new surface
   * @param lump             The lump into which the surface is to be inserted.
   * @return If an error occured, the return value will be NULL.  If the lump
   *         was split, the return value will be the new lump.  Otherwise it
   *         will be the passed lump.
   */
  Lump* insert_surface( Surface* surf_to_copy, Lump* lump );
  
  /** Partition a lump.
   *
   * Create a new surface in the lump from the passed facet patch.
   * This function may partition the lump if the surface divides
   * a shell of the lump into two distinct regions.
   * @param facets  The geometry from which to create the new surface
   * @param lump    The lump into which the surface is to be inserted.
   * @return If an error occured, the return value will be NULL.  If the lump
   *         was split, the return value will be the new lump.  Otherwise it
   *         will be the passed lump.
   */
  Surface* insert_surface( DLIList<CubitFacet*>& facets, Lump* lump );
  
  
  /** Undo surface partitioning
   *
   * Remove a previously inserted split curve.
   * @param curve The curve to remove.
   * @param dead_surfs If non-null, will be populated with pointers to
   *                   any surfaces destroyed by this operation.
   * @return The surviving surface.
   */
  Surface* remove_curve( PartitionCurve* curve,
                         PartitionSurface* dead_surfs[2] = 0 );			
  
  /** Undo lump partitioning
   *
   * Remove a previously inserted split surface.
   * @param surface The surface to remove.
   * @return The surviving lump.
   */
  Lump* remove_surface( PartitionSurface* surface );

  /** Return layer in the TopologyBridge graph that this engine handles */
  int layer() { return SUBCOMP_PARTITION_LAYER; }
  
  //static void destroy_facet( CubitFacetData* facet );
  
  /** Add SubEntitySet to ID map (used for save/restore */
  CubitStatus add_to_id_map( SubEntitySet* set, int unique_id );

  /** Remove SubEntitySet from ID map (used for save/restore */
  CubitStatus remove_from_id_map( SubEntitySet* set, int unique_id );

  /** Retreive SubEntitySet from ID map given ID. */
  SubEntitySet* get_from_id_map( int unique_id );

  /** Retreive a PartitionEntity given ID pair.
   *
   * PartitionEntities are identified by the ID of their
   * SubEntitySet and the ID of the entity in that set. 
   * Retreive the PartitionEntity given these IDs.
   *
   * @param set_id The UniqueID of the SubEntitySet.
   * @param geom_id The ID of the PartitionEntity within the
   *                SubEntitySet
   * @param default_set The SubEntitySet to retreive the 
   *                entity from if the set_id == 0.
   * @return The requested PartitionEntit, or NULL if it does
   *         not exist.
   */
  PartitionEntity* entity_from_id( int set_id, int geom_id, 
                                   SubEntitySet& default_set );
  
  
  /** Destroy a surface and child topology.  
   *
   * Destroy a surface with no parent topology, destroying any
   * child topology that is used only by the passed surface.
   */
  CubitStatus destroy_surface( PartitionSurface* surface );
  
  /** Destroy a lump and child topology.  
   *
   * Destroy the passed lump and any child topology used
   * only by the passed lump.
   */
  CubitStatus destroy_lump( PartitionLump* lump );

  static void delete_facet( CubitFacet* );

  CubitStatus delete_solid_model_entities( PartitionBody* body,
                                           BodySM*& real_body );
  CubitStatus delete_solid_model_entities( PartitionSurface* surface,
                                           Surface*& real_surface );
  CubitStatus delete_solid_model_entities( PartitionCurve* curve,
                                           Curve*& real_curve );

  static void curves_from_surfaces( const DLIList<PartitionSurface*>& input_surfs,
                                    DLIList<PartitionCurve*>& output_curves );

  void clean_out_deactivated_geometry();
  
  CubitStatus translate( PartitionBody* ent, const CubitVector& delta );
  CubitStatus rotate( PartitionBody* ent, const CubitVector& axis, double degrees );
  CubitStatus scale( PartitionBody* ent, const CubitVector& factors );
  CubitStatus reflect( PartitionBody* ent, const CubitVector& axis );
  CubitStatus restore_transform( PartitionBody* ent );
  
  CubitStatus translate( PartitionEntity* ent, const CubitVector& delta );
  CubitStatus rotate( PartitionEntity* ent, const CubitVector& axis, double degrees );
  CubitStatus scale( PartitionEntity* ent, const CubitVector& factors );
  CubitStatus reflect( PartitionEntity* ent, const CubitVector& axis );

  CubitStatus notify_transform( TopologyBridge* bridge,
                                const CubitTransformMatrix& xform );
   
  void remove_attributes( DLIList<TopologyBridge*> &bridge_list );
    //remove Composite attributes off of topology bridges
  virtual void remove_modified(DLIList<Surface*> &all_surfs,
    DLIList<Curve*> &all_curves, DLIList<TBPoint*> &all_pts);
    
  void get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                               DLIList<TopologyBridge*> &tbs );
private:

  CubitStatus notify_transform_internal( TopologyBridge* bridge,
                                const CubitTransformMatrix& xform );

  PartitionEngine();

  friend class PartitionLumpImprint;

  /** Find location to insert a curve into a loop.
   *
   * Given a curve and an end point on the curve, find
   * the two coedges in the loop between which the curve
   * should be inserted.
   *
   *@param surface The surface the curve is to be inserted into.
   *@param curve   The curve to be inserted.
   *@param point   The end point of the curve to find the loop intersection for
   *@param previous A CoEdge ending at the passed point (output).
   *@param next     A CoEdge beginning at the passed point (output).
   */
  CubitStatus find_coedges( PartitionSurface* surface,
                            PartitionCurve* curve,
                            PartitionPoint* point,
                            PartitionCoEdge*& previous,
                            PartitionCoEdge*& next );
  
  /** Helper function for find_coedges(..).
   *
   * Use facet connectivity to find next/prev curve in loop.
   *
   *@param surface The owning surface
   *@param edge    A facet edge owned by the surface
   *@param point   The facet point to search around
   *@param backwards The search direction (false for clockwise)
   *@return The next curve in the specified order around the passed facet point.
   */
  PartitionCurve* next_curve_around_point( PartitionSurface *const surface,
                                           CubitFacetEdgeData *const edge,
                                           CubitPointData *const point,
                                           const bool backwards );
  

  /** Remove a point-curve */
  CubitStatus remove_point_curve( PartitionPoint* point );
  
  /** Create a point curve 
   * 
   * Create point-curve topology.
   * After the surface facetting has been updated with the 
   * location of a point-curve, create the point curve given
   * the facet point.
   *
   * Called by insert_point_curve(ParttionSurface*, const CubitVector& )
   */
  PartitionPoint* insert_point_curve( PartitionSurface* surface,
                                      CubitPointData* point );
  CubitStatus insert_point_curve( PartitionSurface* surface,
                                  PartPTCurve* curve,
                                  bool update_topology = true );
                                  
  
  /** Partition a curve using the passed point.
   *
   * Split a curve.
   *
   * Used by restore_from_attrib(), insert_point(PartitionCurve*,double),
   * and insert_point(CubitPointData*)
   *
   * @param curve  The curve to split.
   * @param point  The new point splitting the curve.
   * @return       The new curve.
   */
  PartitionCurve* insert_point( PartitionCurve* curve, PartitionPoint* point );
  
  /** Partition surfaces
   * 
   * Split surfaces by projecting the passed polyline onto the patch
   * of surfaces.
   * @param surfaces     The list of surfaces to split.
   * @param polyline     The curve to imprint the surfaces with.
   * @param new_surfaces Populated with new surfaces created during partitioning.
   * @param polyline_curves Popolated with new curves.
   */
  CubitStatus insert_curve( DLIList<PartitionSurface*>& surfaces,
                            DLIList<CubitVector*>& polyline,
                            DLIList<PartitionSurface*>& new_surfaces,
                            DLIList<PartitionCurve*>& polyline_curves,
                            const double *tolerance_length = NULL);

  /** Insert a new curve into surface topology.
   *
   * Insert a new curve into the topology of the passed surface,
   * splitting the surface if the curve splits a loop in the surface.
   *
   * The surface facet connectivity is used to determine which
   * curves belong with which surface after the split.  This function
   * assumes that the surface facetting has already been updated with
   * the projection of the new curve, and that the resulting facet 
   * edges are associated with the passed curve.
   * 
   * Used by insert_curve(DLIList<CubitFacetEdgeData*>&), 
   * and restore_from_attrib(Surface*).
   *
   * @param surface  The surface to split
   * @param curve    The split curve
   * @return NULL if an error occured, the new surface if the surface
   *         was split, otherwise the passed surface.
   */
  PartitionSurface* insert_curve( PartitionSurface* surface, 
                                  SegmentedCurve* curve );
  
  /** Insert a curve into surface topology.
   *
   * Given the imprint if a curve on the surface facetting, create
   * a new curve and insert it in the surface topology, splitting
   * the surface if necessary.  The surface is determined from the
   * owners of the facets adjacent to the passed edges.  The new
   * curve is returned.
   *
   * If the list of curve segments is not specified, the facet edges
   * will be used to define the curve geometry.  The list of curve
   * segments may optionally be specified to provide a more precise
   * definition of the curve geometry.
   *
   * Used by insert_Curve( DLIList<Surface*>&, DLIList<CubitVector*&, ... );
   */
  SegmentedCurve* insert_curve( DLIList<CubitFacetEdgeData*>& facet_edges,
                                DLIList<CubitVector*>* curve_segments = 0 );
  
  /** Project a point into the surface facetting
   *
   * Project a position onto the facetting of a surface, updating
   * the surface facetting to contain the resulting point.
   *
   * Used by insert_point_curve(Surface*, const CubitVector&), and
   * restore_from_attrib(Surface*)
   *
   * @param surface   The surface to project onto.
   * @param position  The position to project.
   * @return The new facet point added to the surface facetting.
   */
  CubitPointData* project_to_surface( PartitionSurface* surface,
                                      const CubitVector& position );
  
  /** Project a polyline onto a facet patch.
    * 
    * Project a polyline onto a patch of facets.  Modifies facets
    * and updates owners of facets.  The passed facet patch will
    * no longer be valid after this operation.  If necessary, the
    * facets should be re-aquired from the owning geometry, which
    * will have been updated appropriately.
    *
    * Used by insert_curve(DLIList<Surface*>&, DLIList<CubitVector*>&, ...),
    * and restore_from_attrib(Surface*)
    *
    * @param facets       The facet patch to project onto.
    * @param polyline_in  The polyline to project.
    * @param polyline_out The resulting polyline created in the
    *                     facetting.
    * @param polyline_pts The points in the facetting corresponding to
    *                     each input position in polyline_in.  List
    *                     may contain NULL values.
    */
  CubitStatus project_to_surface( DLIList<CubitFacetData*>& facets,
                                  DLIList<CubitVector*>& polyline_in,
                                  DLIList<CubitFacetEdgeData*>& polyline_out,
                                  DLIList<CubitPointData*>& polyline_pts,
                                  const double *tolerance_length = NULL);
  
  /** Insert a surface into a lump
   *
   * Update lump topology for the creation of a new surface.  Splits
   * the lump if necessary.  
   *
   * Used by insert_surface(Surface*,Lump*) and restore_from_attrib(Lump*)
   *
   * @param lump  The lump to insert the surface in
   * @param surf  The surface to insert
   * @return      NULL if an error occured, the new lump if the passed
   *              lump was split, or the passed lump otherwise.
   */
  PartitionLump* insert_surface( PartitionLump* lump, PartitionSurface* surf );
  
  
  /** Find or create the owning body for a new PartitionPoint */
  PartitionBody* make_body( PartitionPoint* pt );
  /** Find or create the owning body for a new PartitionCurve */
  PartitionBody* make_body( PartitionCurve* curve );
  /** Find or create the owning body for a new PartitionSurface */
  PartitionBody* make_body( PartitionSurface* surf );
  /** Find or create the owning body for a new PartitionLump */
  PartitionBody* make_body( PartitionLump* lump );
  /** Common code used by make_body(..) methods. */
  PartitionBody* make_body_internal( TopologyBridge* child_ptr );

  /** Replicate topology in partition layer */
  PartitionLump* replace_lump( Lump* lump );
  /** Remove partition layer topolgy and restore real topology */
  Lump* restore_lump( PartitionLump* lump );
  /** Helper for destroying PartitionLump topology*/
  CubitStatus destroy_shell( PartitionShell* shell );


public: /* temporarily public for debugging */
  /** Replicate topology in partition layer */
  SubSurface* replace_surface( Surface* surface );
private:
  /** Remove partition layer topolgy and restore real topology */
  Surface* restore_surface( SubSurface* surface );

  /** Replicate topology in partition layer */
  SubCurve* replace_curve( Curve* curve );
  /** Remove partition layer topolgy and restore real topology */
  Curve* restore_curve( SubCurve* curve );
  
  /** Replicate topology in partition layer */
  PartitionPoint* replace_point( TBPoint* point );
  /** Remove partition layer topolgy and restore real topology */
  TBPoint* restore_point( PartitionPoint* point );
  
  /** Split a surface
   * 
   * Split surface and move loops to new surface.
   * When inserting a curve results in a loop being split 
   * (or a new hole-loop being created) this method is called.
   * This method uses facet connectivity and facet-owenr 
   * associativity to move the appropriate surface facets
   * and PartitionLoops to the new surface.
   *
   * Used by insert_curve(PartitionSurface*, SegmentedCurve*).
   *
   *@param surface   The surface to split.
   *@param new_curve A coedge owning the curve that caused the split.
   *@return NULL if an error occured, the new surface if one was
   *        created, or the passed surface if the surface isn't to
   *        be split yet.
   */
  PartitionSurface* split_surface( PartitionSurface* surface, 
                                   PartitionCoEdge* new_curve );
  
  /** Split a lump
   *
   * Split a lump and move shells to new lump.
   * If inserting a surface could potentially result in the
   * lump being partitioned, this method is called to attempt
   * to split the shell and move the appropriate shells to the
   * new lump.
   *
   * Called by insert_surface( PartitionLump*, PartitionSurface*)
   *
   *@param shell  The shell that is potentially to be split.
   *@return NULL if an error occured, a pointer to the new lump
   *        if the shell was split.  The orignal shells owner otherwise.
   */
  PartitionLump* split_lump( PartitionShell* shell );
  
  /** Split a shell
   *
   * Check if a shell can be split, and if so split it.
   *
   * Called by split_lump(PartitionShell*) and remove_surface(PartitionSurface*)
   */
  PartitionShell* split_shell( PartitionShell* shell );
  
  /** Put non-manifold surfaces in appropriate shell.
   *
   * After a shell is split, figure out which of the new shells
   * each surface in the passed list of non-manifold surfaces belongs 
   * in, and add it to that shell.
   *
   * Called by split_shell(PartitionShell*)
   *
   * @param surfaces  A list of non-manifold surfaces with co-surfaces
   *                  already created, but not beloning to any shell.
   * @param shell1    One of the two shells.
   * @param shell2    The other of the two shells.
   */
  void insert_nonmanifold_surfaces( DLIList<PartitionSurface*>& surfaces,
                             PartitionShell* shell1, PartitionShell* shell2);
  
  /** Test if a surface is on the inside of a shell.
   *
   * Given a surface which shares at least one curve with the
   * boundary of the shell (i.e. two non-manifold surfaces in the
   * shell), determine if the surface is inside or outside of the
   * shell.
   *
   * This is a helper function for insert_nonmanifold_surfaces(..)
   */
  CubitStatus inside_shell( PartitionShell* const shell, 
                            PartitionSurface* const surf, 
                            bool& result );
  
  /** Save partition curve geometry as attributes.*/
  CubitStatus save_curves( SubEntitySet* curve_set );
  /** Save surface partitions as attributes. */
  CubitStatus save_surfaces( SubEntitySet* surface_set );
  /** Save lump partitions as attributes. */
  CubitStatus save_lumps( SubEntitySet* lump_set );
  
  /** Restore curve partitions from attributes */
  CubitStatus restore_from_attrib( Curve* partitioned_curve );
  /** Restore surface partitions from attributes */
  CubitStatus restore_from_attrib( Surface* partitioned_surface );
  /** Restore lump partitions from attributes */
  CubitStatus restore_from_attrib( Lump* partitioned_lump );
  
  static PartitionEngine* instance_;

  /** Map for getting a SubEntitySet given its unique ID */  
  std::map<int,SubEntitySet*> uniqueIdMap;

  DLIList<TopologyBridge*> deactivatedList;
};

#endif
