/** \class ChollaEngine
//-------------------------------------------------------------------------
// Filename      : ChollaEngine.hpp 
//
// Purpose       : For a list of facets, builds a topology of
//                 surfaces, curves and vertices 
//
// Special Notes : 
//
// Creator       : Steven J. Owen
//
// Date          : 4/29/2001
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------
*/
#ifndef CHOLLAENGINE_HPP
#define CHOLLAENGINE_HPP

#include "DLIList.hpp"

class FacetEntity;
class ChollaCurve;
class ChollaSurface;
class ChollaPoint;
class CubitPoint;
class CubitFacetEdge;
class CubitFacet;
class GMem;

class ChollaEngine
{
public:
  //! Constructor
  ChollaEngine( );
  //! Constructor
  ChollaEngine(DLIList<FacetEntity*> &face_list,
               DLIList<FacetEntity*> &edge_list,
               DLIList<FacetEntity*> &point_list );
  //! Constructor
  ChollaEngine(DLIList<CubitFacet*> &facet_list,
               DLIList<CubitFacetEdge*> &edge_list,
               DLIList<CubitPoint*> &point_list );
  
    //! Destructor
  ~ChollaEngine();
  void delete_me();

  //! Creates the geometry infrastructure based on the given mesh data.
  CubitStatus create_geometry(CubitBoolean use_feature_angle = CUBIT_TRUE,
                              double angle = 135.0,
                              int interp_order = 0,
                              CubitBoolean smooth_non_manifold = CUBIT_TRUE,
                              CubitBoolean split_surfaces = CUBIT_FALSE);

  void get_surfaces( DLIList<ChollaSurface *> & cholla_surface_list )
    { cholla_surface_list += chollaSurfaceList; }
  void get_curves( DLIList<ChollaCurve *> & cholla_curve_list )
    { cholla_curve_list += chollaCurveList; }
  void get_points( DLIList<ChollaPoint *> & cholla_point_list )
    { cholla_point_list += chollaPointList; }

  void delete_eval_tools();
    //! delete the eval tools associated with the cholla geom entities
  void delete_eval_tools_but_not_facets();
    

    //! Create the new features by cracking the facets at feature edges.
    //! Creates new edges and points and updates connectivity

  static CubitStatus make_features( DLIList<CubitFacetEdge *> &feature_edge_list,
                             CubitBoolean split_surfaces );

#ifdef ALPHA_CABLE
  CubitStatus split_surface_at_edges( ChollaSurface* owning_surface,
                                      DLIList<CubitFacetEdge *> &feature_edge_list,
                                      DLIList<ChollaCurve*> &cholla_curves,
                                      DLIList<ChollaSurface*> &cholla_surfaces);
    // given a list of edges, crack the surface.
#endif
    //! assumes that valid CubitPoint* are supplied
    //! marks the CubitPoints as features to manually break up curves
    // added by Ved Vyas 8-2003
  static void mark_features (DLIList<CubitPoint*> &feature_points);

    //! fix the control points so they are C-zero          
  static CubitStatus fix_geometry( CubitBoolean smooth_non_manifold,
                                   CubitBoolean split_surfaces,
                                   CubitBoolean use_feature_angle,
                                   double min_dot,
                                   DLIList <CubitFacetEdge *> &feature_edge_list );

    //! ensure facets are all oriented consistently
    //Cubit is not using this anymore.
  static CubitStatus check_all_facet_orientations( DLIList <CubitFacet*> &facet_list,
                                                  CubitBoolean do_flip = CUBIT_FALSE);

    //! determine orienation of curve w.r.t. the surface    
  static CubitStatus determine_curve_orientation( ChollaSurface *chsurf_ptr,
                                                  ChollaCurve *chcurv_ptr,
                                                  CubitSense & orientation );

  static CubitStatus get_facets(GMem& gmem, DLIList<CubitFacet*> &facet_list, 
                                DLIList<CubitPoint*> &dl_point_list);

    //set the ChollaEngine to actually flip facets, or to just set a flag
    // when facets need to be reoriented
  void set_flip_flag(CubitBoolean flip)
    {
      doFlip = flip;
    }
    
  void print_me();
  void dump( char *filename, double angle );
   
private:

  DLIList<FacetEntity*> faceList;
  DLIList<FacetEntity*> edgeList;
  DLIList<FacetEntity*> pointList;
  DLIList<ChollaSurface*> chollaSurfaceList;
  DLIList<ChollaCurve*> chollaCurveList;
  DLIList<ChollaPoint*> chollaPointList;

  DLIList<ChollaCurve*> *hashCurveArray;
  int hashCurveSize;
  DLIList<ChollaPoint*> *hashPointArray;
  int hashPointSize;
    //boolean to determine whether flip the facets (when needed) or set
    // the isBackwards flag.
  CubitBoolean doFlip;

  void set_up_tool_datas( );
  void delete_tool_datas( );
  
    //! returns the dimension of the facet entity.
  int facet_dimension(FacetEntity *facet_ptr);
   
    //! creates the correct facetsurfaces for all the element blocks.
  CubitStatus create_volume_boundaries( DLIList<ChollaSurface*> &facet_surface_sheets,
                                        CubitBoolean use_feature_angle,
                                        double min_dot,
                                        CubitBoolean split_surfaces );

    //! creates the correct blockcurves for all the surfaces  
  CubitStatus create_surface_boundaries( DLIList<ChollaSurface*> &cholla_surface_list,
                                         DLIList<ChollaCurve*> &cholla_curve_list,
                                         CubitBoolean use_feature_angle,
                                         double min_dot );
  
    //! creates the correct blockpoints for all the curves
  CubitStatus create_curve_boundaries( DLIList<ChollaCurve*> &cholla_curve_list,
                                       DLIList<ChollaPoint*> &cholla_point_list );
  
    //! sorts a edge into its correct curve based on its associated
    //! blocks and sidesets.  Creates a new facet curve if necessary.
  CubitStatus classify_edge( FacetEntity *edge_ptr,
                             DLIList<ChollaCurve*> &cholla_curve_list,
                             ChollaSurface *fsm_ptr );
    
    //! sorts a node into correct point based on its associated
    //! curve.  Creates a new facet point if necessary
  CubitStatus classify_point(CubitPoint *point_ptr,
                             DLIList<ChollaPoint*> &cholla_point_list,
                             ChollaCurve *fcm_ptr ); 
       
    //! build the CUBIT geometry based on the Facet entity class lists
  CubitStatus build_eval_tools( DLIList<ChollaSurface*> &cholla_surface_list,
                                DLIList<ChollaCurve*> &cholla_curve_list,
                                int interp_order,
                                CubitBoolean use_feature_angle,
                                double min_dot,
                                CubitBoolean smooth_non_manifold ,
                                CubitBoolean split_surfaces );
    //! build the CUBIT geometry based on the Facet entity class lists
    //! From the facet curve list, create geometric curves for each
  CubitStatus build_curve_eval_tools( DLIList<ChollaCurve*> &cholla_curve_list,
                                    int interp_order );

    //! From the facet surface list, create geometric surface,
    //! loops and coedges for each surface in the list
  CubitStatus build_surface_eval_tools( DLIList<ChollaSurface*> &cholla_surface_list,
                                      int interp_order, double min_dot );

    //! functions for hashing curves - to speed up edge classification
  CubitStatus init_hash_curves( );
  void delete_hash_curves( );
  int get_curve_hash_key( DLIList<ChollaSurface*> *bsm_list_ptr );

    //! functions for hashing points - to speed up node classification
  CubitStatus init_hash_points( );
  void delete_hash_points( );
  int get_point_hash_key( DLIList<ChollaCurve*> *bsm_list_ptr );

    //! fix the edge control points and the normals so they are conforming
    //! (or non-conforming) accross curves
  CubitStatus clean_geometry( CubitBoolean smooth_non_manifold,
                            CubitBoolean split_surfaces,
                            CubitBoolean use_feature_angle,
                            double min_dot,
                            DLIList <ChollaCurve *> &cholla_curve_list );


    //! make sure the facets are oriented consistently
  static CubitStatus check_facet_orientation( CubitFacet *facet, 
                                              CubitBoolean do_flip, 
                                              int &nfacets, int mydebug = 0 );

    //! split the surface facets at a point if required based on features.  
    //! Create new points and edges as needed
  static CubitStatus crack_open_point( CubitPoint *point_ptr, int mydebug );

    //! does the same as crack_open_point, but does not create new facet entities
    //! at surface boundaries.  Instead, it creates facet boundary tooldatas
    //! to hold the additional control point and normal data
  static CubitStatus insert_discontinuity_at_point( CubitPoint *point_ptr );

    //! recursive function - returns a list of facets adjacent to
    //! a point that are bounded by feature edges
  static CubitStatus get_facets_at_point( CubitPoint *point_ptr,
                                   CubitFacet *facet_ptr,
                                   DLIList<CubitFacet *> &facet_list,
                                   DLIList<CubitFacetEdge *> &feature_edge_list );
    //! fix control points on a single edge and its partners
  static CubitStatus fix_edge( CubitFacetEdge *edge_ptr, 
                        DLIList<CubitFacet *> &update_facet_list );
  static CubitStatus fix_split_edge( CubitFacetEdge *edge_ptr, 
                        DLIList<CubitFacet *> &update_facet_list );

    //! fix normals at non-manifold edges to maintain
    //! continuity across facets that don't meet feature angle criteria
  static CubitStatus fix_split_non_manifold_edge( CubitFacetEdge *edge_ptr,
                                     double min_dot,
                                     DLIList <CubitPoint *> &changed_points );
  static CubitStatus fix_non_manifold_edge( CubitFacetEdge *edge_ptr,
                                     double min_dot,
                                     DLIList <CubitPoint *> &changed_points );

    //! update the edge control points at all edges connected to points
    //! in point_list
  static CubitStatus update_edges_at_points( CubitBoolean split_surfaces,
                                      DLIList<CubitPoint *> &point_list,
                                      DLIList<CubitFacet *> &facet_update_list,
                                      double mindot );
    //! set the normals on two points so they are the same
  static CubitStatus merge_normals( CubitPoint *pt0, CubitPoint *pt1);

};

#endif


