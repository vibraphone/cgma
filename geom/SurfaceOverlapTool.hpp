//- Class:       SurfaceOverlapTool
//- Description: Utilities to debug imprinting/merging problems
//- Owner:       Steve Storm
//- Created:     22 October 1999

#ifndef SurfaceOverlapTool_HPP
#define SurfaceOverlapTool_HPP

#include <map>
#include "GMem.hpp"
#include "AnalyticGeometryTool.hpp"
#include "CubitGeomConfigure.h"
#include "CurveOverlapFacet.hpp"
#include "SurfaceOverlapFacet.hpp"
#include "AbstractTree.hpp"

template <class X> class DLIList;
class RefFace;

//! Tool to find surfaces and curves that overlap.
class CUBIT_GEOM_EXPORT SurfaceOverlapTool
{

public :   

   //! From the specified surfaces, find surfaces that overlap. 
   CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list, 
                                          DLIList<RefEntity*> &faces_to_draw,
                                          bool filter_slivers = false );

   //! From the specified bodies, find surfaces that overlap. 
   CubitStatus find_overlapping_surfaces( DLIList<Body*> &body_list,
                                          DLIList<RefEntity*> &faces_to_draw,
                                          bool filter_slivers = false );

   //! From the specified surfaces, find surfaces that overlap. 
   CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                          DLIList<RefFace*> &ref_face_list1,
                                          DLIList<RefFace*> &ref_face_list2,
                                          DLIList<RefEntity*> &faces_to_draw,
                                          CubitBoolean show_messages = CUBIT_FALSE,
                                          bool filter_slivers = false);

   //! From the specified bodies, find surfaces that overlap. 
   CubitStatus find_overlapping_surfaces( DLIList<Body*> &body_list,
                                          DLIList<RefFace*> &ref_face_list1,
                                          DLIList<RefFace*> &ref_face_list2,
                                          DLIList<RefEntity*> &faces_to_draw,
                                          CubitBoolean show_messages = CUBIT_FALSE,
                                          bool filter_slivers = false);

   //! From the specified bodies, find surfaces that overlap. 
   CubitStatus find_candidate_surfaces_for_imprinting( DLIList<BodySM*> &body_list,
                                          DLIList<Surface*> &surface_list1,
                                          DLIList<Surface*> &surface_list2,
                                          double overlap_tol = -1.0,
                                          bool filter_slivers = false );

   // Searches for surfaces that overlap each other and are good
   // candidates for imprinting.  It can even find those with gaps
   // between them (according to settings).
   // First method is for interactive use.  They print messages and
   // draw the surfaces as it finds them, as dictated by settings.
   // The second method is for non-interactive use.  It simply returns the
   // surface pairs in lists.  Note that the same surface can occur 
   // multiple times in the lists.  The first list will contain sequential
   // surfaces - ie., if surface 1 overlaps with 2 other surfaces, list1
   // will contain surface 1, then surface 1, then the next surface, 
   // etc..  If show_messages is CUBIT_TRUE, you will see "Finding surface
   // overlap...", and groupResults, listPairs and displayPairs will be
   // used, otherwise they will be false (default).

   //! From the specified surfaces, find overlapping curves.
   CubitStatus find_overlapping_curves( DLIList<Surface*> &surf_list,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_point_imprint_map);

   //! From the specified curves, find those that overlap. 
   CubitStatus find_overlapping_curves( DLIList<Curve*> &curve_list,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_point_imprint_map, double overlap_tol=-1.0);

   //! From the specified bodies, find curves that overlap. 
   CubitStatus find_overlapping_curves( DLIList<Body*> &bodies,
                   std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                   double maxgap = -1);

   //! From the specified bodies, find curves that overlap. 
   CubitStatus find_overlapping_curves( DLIList<BodySM*> &bodies,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_vertex_imprint_map,
                                double overlap_tol);

   //! From the specified surfaces, find curves that overlap. 
   CubitStatus find_overlapping_curves( DLIList<RefFace*> &faces,
				        std::multimap<RefEdge*, 
                        RefEdge*> &overlapping_edge_map,
                        double maxgap = -1);


   //! From the specified curves, find curves that overlap. 
   CubitStatus find_overlapping_curves( DLIList<RefEdge*> &edgelist,
				        std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                        double maxgap = -1);

   //! From the specified curves, find curves that overlap. 
   CubitStatus find_overlapping_curves( DLIList<RefEdge*> &edges1, DLIList<RefEdge*> &edges2,
					std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                    double maxgap = -1);

  //! Checks to if two facets overlap
  CubitBoolean check_overlap( SurfaceOverlapFacet *facet1, SurfaceOverlapFacet *facet2, const double facet_compare_tol );

   //! Checks to see if two curves overlap.  Reuses graphic facets.
   CubitBoolean check_overlap( RefEdge *edge1, RefEdge *edge2, 
               std::map<RefEdge*, DLIList<CurveOverlapFacet*>* > *facet_map, 
               double *overlap_tol = NULL );

   //! Checks to see if two surfaces overlap.
   CubitBoolean check_overlap( RefFace *ref_face_ptr1,
                               RefFace *ref_face_ptr2,
                               CubitBoolean abort,
                               CubitBoolean draw_overlap = CUBIT_FALSE,
                               double *overlap_area = NULL );

   //! Checks to see if two curves overlap.  Reuses graphic facets.
   CubitBoolean check_overlap( Curve *curve1, Curve *curve2, 
               std::map<Curve*, DLIList<CurveOverlapFacet*>* > *facet_map, 
               std::multimap<BodySM*, CubitVector > *body_point_imprint_map = NULL,
               double overlap_tol = -1.0);

   //! Checks for boundary contact between surfaces
   CubitBoolean check_boundary_contact(  
                                DLIList<SurfaceOverlapFacet*> *facet_list1, 
                                DLIList<SurfaceOverlapFacet*> *facet_list2, 
                                AbstractTree<SurfaceOverlapFacet*> *a_tree, 
                                const double facet_compare_tol, 
                                const double tmp_overlap_tol );

   
  // swap surfaces and facets to keep the smaller surface first
  CubitBoolean check_size_and_swap_surfs( 
                                Surface *&tmp_surf1, 
                                Surface *&tmp_surf2, 
                                const double tolerance,
                                const double facet_tol,
                                DLIList<SurfaceOverlapFacet*> *&facet_list1, 
                                DLIList<SurfaceOverlapFacet*> *&facet_list2, 
                                std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* > *&a_tree_map,
                                AbstractTree<SurfaceOverlapFacet*> *&a_tree2 );

  // extract facets if not already available 
   CubitBoolean extract_surf_facets( 
                               Surface *surface1, 
                               Surface *surface2, 
                               std::map<Surface*, DLIList<SurfaceOverlapFacet*>* > *facet_map,
                               const double tolerance,
                               const double facet_tol,
                               DLIList<SurfaceOverlapFacet*> *&facet_list1,
                               DLIList<SurfaceOverlapFacet*> *&facet_list2 );

   //! Determine tolerance for checking intersection between surfaces
   CubitBoolean calculate_tolerances_for_surf_intersection( 
                                Surface *tmp_surf1, 
                                Surface *tmp_surf2, 
                                DLIList<SurfaceOverlapFacet*> *facet_list1, 
                                DLIList<SurfaceOverlapFacet*> *facet_list2, 
                                std::map<Surface*, double > *area_map, 
                                double &facet_compare_tol, 
                                double &overlap_tolerance );

   //! Checks for overlap betwen facet_lists of two surfaces
   CubitBoolean check_overlap(  DLIList<SurfaceOverlapFacet*> *facet_list1, 
                                DLIList<SurfaceOverlapFacet*> *facet_list2, 
                                AbstractTree<SurfaceOverlapFacet*> *a_tree, 
                                const double facet_tol, 
                                const double tolerance );   

   //! Checks to see the two groups of facets overlap. 
   CubitBoolean check_overlap( DLIList<SurfaceOverlapFacet*> *facet_list1,
                               DLIList<SurfaceOverlapFacet*> *facet_list2,
                               AbstractTree<SurfaceOverlapFacet*> *a_tree, 
                               CubitBoolean abort, 
                               CubitBoolean draw_overlap,
                               double *overlap_area );

   //! Draws the overlapping surface pair.
   CubitBoolean draw_overlapping_surface_pair( RefFace *ref_face_1,
                                               RefFace *ref_face_2);


   //! List the settings used for surface overlap detection
   void list_settings();

   unsigned short get_facet_ang_tol();
   void set_facet_ang_tol( unsigned short angle_tol );


   //Initialize all settings in this class
   static void initialize_settings();

   static int get_facet_ang_tol_setting();
   static void set_facet_ang_tol_setting( int angle_tol );
   
   static double get_facet_abs_tol();
   static void set_facet_abs_tol( double absolute_tol );

   static double get_gap_min();
   static void set_gap_min( double val );

   static double get_gap_max();
   static void set_gap_max( double val );

   static double get_angle_min();
   static void set_angle_min( double val );

   static double get_angle_max();
   static void set_angle_max( double val );

   static double get_overlap_tolerance();
   static void set_overlap_tolerance( double val );

   static int get_normal_type(); //
   static void set_normal_type( int type ); //

   static CubitString get_normal_type_setting();
   static void set_normal_type_setting( CubitString type );

   static CubitBoolean get_group_results(); //
   static void set_group_results( CubitBoolean setting ); //

   static CubitBoolean get_list_pairs(); //
   static void set_list_pairs( CubitBoolean setting ); //

   static CubitBoolean get_display_pairs(); //
   static void set_display_pairs( CubitBoolean setting ); //

   static int get_group_results_setting();
   static void set_group_results_setting( int setting );

   static int get_list_pairs_setting();
   static  void set_list_pairs_setting( int setting );

   static  int get_display_pairs_setting();
   static  void set_display_pairs_setting( int setting );

   static CubitBoolean get_imprint();
   static void set_imprint( CubitBoolean setting );

   static CubitBoolean get_check_within_bodies();
   static void set_check_within_bodies( CubitBoolean setting );

   static CubitBoolean get_check_across_bodies();
   static void set_check_across_bodies( CubitBoolean setting );
   CubitBoolean get_skip_facing_surfaces();
   void set_skip_facing_surfaces( CubitBoolean setting );
   static SurfaceOverlapTool* instance();
     // Returns a static pointer to unique instance of this class.
   static void delete_instance()
   {
     if(instance_)
       delete instance_;
     instance_ = NULL;
   };
   
   ~SurfaceOverlapTool();
     //- Destructor.

protected :
   
private :



  CubitStatus find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                         DLIList<RefFace*> &ref_face_list1,
                                         DLIList<RefFace*> &ref_face_list2,
                                         DLIList<RefEntity*> &pair_list,
                                         CubitBoolean list_pairs,
                                         int prog_step,
                                         bool filter_slivers = false);
  // Private workhorse function for finding the overlapping surfaces. Here 
  // pair_list is a list of all the surfaces that were found that overlap.
  // If prog_step == -1 it does not step the progress bar. 

  double find_area_overlap( SurfaceOverlapFacet *facet1, SurfaceOverlapFacet *facet2, const double facet_compare_tol );

  
  CubitBoolean check_surfs_for_imprinting( Surface *surface1,
                               Surface *surface2,
              std::map<Surface*, DLIList<SurfaceOverlapFacet*>* > *facet_map,
              std::map<Surface*, double > *area_map,
              std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* > *a_tree_map,
              double overlap_tol=-1.0); 

                                
   //CubitStatus draw_facets( GMem* gMem, int color = 2 );
   //void draw_triangle( Triangle3 &T, int color = 2 );

   CubitStatus imprint(DLIList<RefFace*> &ref_face_list1,
                       DLIList<RefFace*> &ref_face_list2);
   //- Imprints the bodies containing the lists of surfaces.  Also
   //- uses the edges imprinting rather than the normal imprint.

   int num_descendants(Body *body);
    //- Counts the total number of descendants of the body.
    //- Only counts Ref* entities.
  
   CubitStatus copy_edges_in_list(DLIList<RefEdge*> &old_list,
                                  DLIList<RefEdge*> &new_list );
    //- copy's each of the edges in the old list, and puts the
    //- stand alone edges in the new list.
    //- NOTE: THESE ARE NEWLY CREATED EDGES AND SHOULD BE DELETED
    //-       BY THE CALLING FUNCTION.

   CubitStatus delete_edges_in_list(DLIList<RefEdge*> &edge_list);
    //- Deletes the RefEdges in the list via GeometryQueryTool.

   static SurfaceOverlapTool* instance_;
     // Static pointer to unique instance of this class
   
   SurfaceOverlapTool();
     //- Constructor for the SurfaceOverlapTool object

   /*! The angular tolerance indicates the maximum angle between 
       normals of adjacent surface facets. The default angular 
       tolerance is 15° - consider using a value of 5° . This will 
       generate a more accurate facetted representation of the 
       geometry for overlap detection.  */
   static  unsigned short facetAngTol;

  /*!  The distance tolerance means the maximum actual distance 
       between the generated facets and the surface. This value is 
       by default ignored by the facetter - consider specifying a 
       reasonable value here for more accurate results */
   static  double facetAbsTol;

   /*! The overlap algorithm will search for surfaces that are 
   within a distance from the gapMin to gapMax. */
   static  double gapMin;
   static  double gapMax;

   //@{
   /*! Angle comparison for determining if facets are coplanar.  If
   they are within angleMin and angleMax they are considered to be 
   coplanar. */
   static  double angleMin;
   static  double angleMax;
   //@}
  
   //! Searches for overlapping surface that have surface normals
   //! in any, opposite, or the same direction.
   static  int normalType; // 1=any, 2=opposite, 3=same
  
   //! Controls if surfaces are put in a group or not.
   static  CubitBoolean groupResults;
   
   //! Prints out the lists of overlapping surface pairs.
   static  CubitBoolean listPairs;

   //! Draws the lists of overlapping surface pairs.
   static  CubitBoolean displayPairs;
  
   //! Imprint overlapping surfaces parts onto one another. 
   static  CubitBoolean imprintResults;
  
   //! The area threshold that 2 surfaces have to actually overlap.
   //! Prevents detection of sliver overlaps. 
   static  double overlapTolerance;

   //! if the normalType is 2 skip surfaces that face each other
   static CubitBoolean skipFacingSurfaces;

   //! By default this tool will not search for overlapping 
   //! pairs within bodies - only between different bodies. 
   //! Turn this setting on to search for pairs within bodies. 
   //! Note however that this will slow the algorithm down. 
   static  CubitBoolean checkWithinBodies;
  
   //! If true, by default, will try to find surfaces of body A
   //! that overlap surfaces of body B.  Setting this flag to 
   //! false will make the tool only check for surfaces that 
   //! overlap in the same body. 
   static  CubitBoolean checkAcrossBodies;
};

inline
void SurfaceOverlapTool::set_facet_ang_tol( unsigned short val )
{facetAngTol=val;}

inline
double SurfaceOverlapTool::get_facet_abs_tol()
{return facetAbsTol;}
inline
void SurfaceOverlapTool::set_facet_abs_tol( double val )
{facetAbsTol=val;}

inline
double SurfaceOverlapTool::get_gap_min()
{return gapMin;}
inline
void SurfaceOverlapTool::set_gap_min( double val )
{gapMin=val;}

inline
double SurfaceOverlapTool::get_gap_max()
{return gapMax;}
inline
void SurfaceOverlapTool::set_gap_max( double val )
{gapMax=val;}

inline
double SurfaceOverlapTool::get_angle_min()
{return angleMin;}
inline
void SurfaceOverlapTool::set_angle_min( double val )
{angleMin=val;}

inline
double SurfaceOverlapTool::get_angle_max()
{return angleMax;}
inline
void SurfaceOverlapTool::set_angle_max( double val )
{angleMax=val;}

inline
double SurfaceOverlapTool::get_overlap_tolerance()
{return overlapTolerance;}
inline
void SurfaceOverlapTool::set_overlap_tolerance( double val )
{overlapTolerance=val;}

inline
int SurfaceOverlapTool::get_normal_type()
{return normalType;}
inline
void SurfaceOverlapTool::set_normal_type( int type )
{normalType=type;}

inline
CubitBoolean SurfaceOverlapTool::get_group_results()
{return groupResults;}
inline
void SurfaceOverlapTool::set_group_results( CubitBoolean setting )
{groupResults=setting;}

inline
CubitBoolean SurfaceOverlapTool::get_list_pairs()
{return listPairs;}
inline
void SurfaceOverlapTool::set_list_pairs( CubitBoolean setting )
{listPairs=setting;}

inline
CubitBoolean SurfaceOverlapTool::get_display_pairs()
{return displayPairs;}
inline
void SurfaceOverlapTool::set_display_pairs( CubitBoolean setting )
{displayPairs=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_imprint()
{return imprintResults;}
inline
void SurfaceOverlapTool::set_imprint( CubitBoolean setting )
{imprintResults=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_check_within_bodies()
{return checkWithinBodies;}
inline
void SurfaceOverlapTool::set_check_within_bodies( CubitBoolean setting )
{checkWithinBodies=setting;}
inline
CubitBoolean SurfaceOverlapTool::get_check_across_bodies()
{return checkAcrossBodies;}
inline
void SurfaceOverlapTool::set_check_across_bodies( CubitBoolean setting )
{checkAcrossBodies=setting;}  
inline
CubitBoolean SurfaceOverlapTool::get_skip_facing_surfaces()
{return skipFacingSurfaces;}
inline
void SurfaceOverlapTool::set_skip_facing_surfaces( CubitBoolean setting )
{skipFacingSurfaces=setting;}  
#endif



