#ifndef COMPOSITE_TOOL_HPP
#define COMPOSITE_TOOL_HPP

#include "CubitDefines.h"
#include "BasicTopologyEntity.hpp"
#include "GeometryEntity.hpp"

#include "CompositeSurface.hpp"
#include "CompositeCurve.hpp"
#include "CompositePoint.hpp"

class RefEdge;
class RefFace;
class RefEdge;
class RefFace;
template <class X> class DLIList;
class Loop;

class CompositeTool
{
  friend class PartitionTool;
  public:
  
    static CompositeTool* instance();
    
    virtual ~CompositeTool();
    
    virtual RefEdge*   composite( DLIList<RefEdge*>& list,
                                  RefVertex* keep_vertex = NULL,
                                  RefEdge* survivor = NULL );
    virtual RefFace*   composite( DLIList<RefFace*>& list,
                                  RefFace* survivor = NULL );
    virtual RefVolume* composite( DLIList<RefVolume*>& list,
                                  Body** composite_body = 0 );
    virtual Body*      composite( DLIList<Body*>& list );
    //R- the composite, or NULL if error
    //I list
    //I- The list of entities to composite
    //I keep_vertex
    //I- vertex to keep if list forms a closed loop of curves
    //I composite_body
    //I- If non-null, a composite body will be created if necessary,
    //I- and composite_body will be set to point to the new body.
    //- create a composite
    //-
    //- NOTE: The 'survivor' argument is provided to meet
    //- special requirenments in supporting pre-version-8 save
    //- files.  Once support for pre-version-8 save files has
    //- been dropped, this argument will be removed. 
    
    
    
    virtual CubitStatus uncomposite( RefEdge* composite_edge, 
                                     DLIList<RefEdge*>* restored_edges = NULL );
    virtual CubitStatus uncomposite( RefFace* composite_face,
                                     DLIList<RefFace*>* restored_faces = NULL,
                                     bool force_unmerge = false );
    virtual CubitStatus uncomposite( RefVolume* composite_vol,
                                     DLIList<RefVolume*>* restored_vols = NULL );
    virtual CubitStatus uncomposite( Body* composite_body,
                                     DLIList<Body*>* restored_bodies = NULL );
    //Remove a composite entity.
                             
    CubitBoolean isComposite( const BasicTopologyEntity* bte_ptr ) const;
    CubitBoolean isComposite( const GeometryEntity* ge_ptr ) const;
    
    CubitStatus composite( DLIList<RefEdge*>& edges_to_composite, 
                           DLIList<RefEdge*>& new_edge_list,
                           DLIList<RefVertex*>* vertices_to_keep = NULL,
                           double max_angle = 3.15 );
     //R CubitStatus
    //I edges_to_composite
    //I- The list of RefEdges to composite.
    //O new_edge_list
    //O- The list of Composite RefEdges created.
    //I vertices_to_keep
    //I- If this list is passed, no composites will be created which
    //I- would result in the removal of any of these vertices.  Or
    //I- if the curves to be composited form a closed loop, this may
    //I- be used to specify which vertex is kept.
    //I max_angle
    //I- The maximum angle between tangents at the point of intersection
    //I- for which curves can be composited.  The default value (3.15) will
    //I- result in all possible curves being composited.  A value of 0.0
    //I- will result in only C1-continous curves being composited.
    //- Create all possible composites from the passed list of RefEdges,
    //- within the constraints of topological validity, the list of
    //- vertices to keep, and the specified angle.
    
    CubitStatus composite( DLIList<RefFace*>& faces_to_composite,
                           DLIList<RefFace*>& result_faces,
                           double max_angle = 4.0 /*radians*/);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I facets_to_composite
    //I- A list of RefFaces to create composites from.
    //O result_facets
    //O- The resulting composite surfaces.
    //I max_angle
    //I- The largest angle between surfaces to composite over.
    //- Create all possible composites from the passed list of RefFaces,
    //- within the constraints of toopologal validity, and the
    //- specified angle.
    
    
    RefFace* composite( DLIList<RefFace*>& faces_to_composite,
                        DLIList<RefEdge*>& result_edges,
                        DLIList<RefVertex*>* vertices_to_keep = NULL,
                        double max_vertex_angle = 0.26 /*radians*/ );
    //R RefFace*
    //R- The resulting composite, or NULL if an error was encountered.
    //I faces_to_composite
    //I- The list of connected surfaces to composite.
    //O result_edges
    //O- The list of curves composited with the surface.
    //I vertices_to_keep
    //I- The list of vertices that should not be composited over when
    //I- creating composite curves.
    //I max_vertex_angle
    //I- The maximum angle between curves at their common vertex, for
    //I- which the vertex can be composited over.
    //- Create a single composite surface, and as many composite curves
    //- as possible on the boundary of that surface.
    
    CubitStatus composite( DLIList<RefFace*>& faces_to_composite,
                           DLIList<RefFace*>& result_faces,
                            DLIList<RefEdge*>& result_edges,
                           double max_curve_angle = 2.0 /*radians*/,
                           DLIList<RefVertex*>* vertices_to_keep = NULL,
                           double max_vtx_angle = 0.26 /*radians*/ );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I facets_to_composite
    //I- The list of surfaces to create composites from.
    //O result_faces
    //O- The composite surfaces created.
    //O result_edges
    //O- The composite curves created.
    //I max_curve_angle
    //I- The largest angle between surfaces (at a curve) that can
    //I- be composited over.
    //I vertices_to_keep
    //I- The list of vertices that should not be composited over
    //I- when creating composite curves.
    //I max_vtx_angle
    //I- The largest angle between curves (at a vertex) that can
    //I- be composited over.
    //- Create composite surfaces and curves.
    
    
    CubitBoolean okayToComposite( DLIList<BasicTopologyEntity*>& bte_list,
                                  DLIList<BasicTopologyEntity*>* boundary = NULL,
                                  DLIList<BasicTopologyEntity*>* internal = NULL,
                                  bool print_errors = true,
                                  bool force_same_parent = false ) const;
    //This function is public because it is also used in PartitionTool, 
    //and because it is only a series of queries anyway.
    
    static CubitStatus classify_children( DLIList<BasicTopologyEntity*>& bte_list,
                                  DLIList<BasicTopologyEntity*>& boundary_children,
                                  DLIList<BasicTopologyEntity*>& internal_children,
                                  DLIList<BasicTopologyEntity*>& unknown_children );
    
    RefEdge* remove_vertex( RefVertex* vertex, 
                            bool remove_partitions = false,
                            bool update_dag = true,
                            RefEdge* survivor = NULL );
    RefFace* remove_edge( RefEdge* edge, 
                          bool remove_partitions = false,
                          bool update_dag = true,
                          RefFace* survivor = NULL );
    
    RefVolume* remove_face( RefFace* face,
                            bool remove_partitions = false,
                            bool update_dag = true );
    
  protected:
  
    // These functions do nothing as implemented in CompositeTool.
    // They are provided for CompositeToolMesh to overload with 
    // mesh-handling code.
    virtual CubitStatus update_combined_vols( RefVolume* keep, RefVolume* dead );
    virtual CubitStatus update_combined_faces( RefFace* keep, RefEdge* dead, RefFace* delete_face );
    virtual CubitStatus update_combined_edges( RefEdge* keep, RefEdge* dead,int keep_interval,double keep_size,FirmnessType keep_interval_type,SizeIntervalType keep_size_type);
    virtual CubitStatus determine_combined_edges_interval_or_size( RefEdge* edge_1, RefEdge* edge_2,int& result_interval,double& result_size,FirmnessType& interval_keep_type,SizeIntervalType& size_keep_type); 

    CompositeTool();
    
    //CubitStatus make_mergeable( GeometryEntity* bridge1, GeometryEntity* bridge2 );
    
    CubitStatus removeOldTopology( BasicTopologyEntity* composite );
    
    CubitVector tangent( RefEdge* edge_ptr, RefVertex* vtx_ptr ) const;
    //R CubitVector
    //R- The tangent vector
    //I edge_ptr
    //I- The RefEdge to get the tangent of.
    //I vtx_ptr
    //I- The RefVertex to get the tangent at.
    //- Return the tangent of the passed RefEdge at the location of
    //- the passed RefVertex.  The tanget is reversed if the passed
    //- RefVertex is the end vertex of the RefEdge, such that the 
    //- resulting tangent vector points along the edge away from the
    //- vertex.
    

    CubitStatus fast_edge_sort( DLIList<RefEdge*>& edge_list,
                                bool valence2_vertices ) const;
    //R CubitStatus
    //R- CUBIT_SUCCESS if all edges form a single chain.
    //R- CUBIT_FAILURE otherwise
    //I edge_list
    //I- The list of edges to sort
    //O edge_list
    //O- The sorted list of edges
    //- This function sorts edges by common vertices.  The result
    //- will be several consecutive sorted sets of edges, if all
    //- edges do not form a single chain.  
    
    
    CubitStatus restore_merged_point( TBPoint* point, 
                                      DLIList<RefFace*>& update,
                                      bool force = false );
    
    CubitStatus stitch_points( TBPoint* pt1, TBPoint* pt2 );
    CubitStatus stitch_curves( Curve* cv1, Curve* cv2 );
                          
    static CompositeTool* instance_;
    
  private:
  
};

inline CompositeTool* CompositeTool::instance()
{
  return (instance_ != NULL ) ? instance_ : instance_ = new CompositeTool();
}

inline CubitBoolean CompositeTool::
  isComposite( const GeometryEntity* ge_ptr ) const
{
  if( dynamic_cast<const CompositeSurface*>(ge_ptr) ||
      dynamic_cast<const CompositeCurve*>(ge_ptr) ||
      dynamic_cast<const CompositePoint*>(ge_ptr) )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

inline CubitBoolean CompositeTool::
  isComposite( const BasicTopologyEntity* bte_ptr ) const
{
  return bte_ptr ? isComposite(bte_ptr->get_geometry_entity_ptr()) : CUBIT_FALSE;
}


#endif
