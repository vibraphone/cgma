#ifndef FACETBOOLINTERFACE_HPP
#define FACETBOOLINTERFACE_HPP

class OCCBody;
class BodySM;
class OCCSurface;
class OCCCurve;
class Curve;
class CubitBox;

#include <vector>
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "FBStructs.hpp"

class FacetboolInterface {

public:
  
  FacetboolInterface();
  ~FacetboolInterface();

CubitStatus webcut_FB(BodySM *bodysm_ptr,
                                         STD(vector)<double>& cutter_verts,
                                         STD(vector)<int>& cutter_connections,
                                         bool cutter_is_plane,
                                         CubitBoolean delete_bodies,
                                         CubitBoolean &intersects,
                                         DLIList<BodySM*>& results_list);
                                         
  CubitStatus dofacetboolean(DLIList<BodySM*>& body_list, 
                             BodySM*& newBody,
                             bool keep_old,
                             const CubitFacetboolOp op);

  CubitStatus dofacetboolean_subtract(BodySM*& tool_body, 
                                      DLIList<BodySM*>& from_bodies,
                                      DLIList<BodySM*>& new_bodies,
                                      bool keep_old,
                                      bool* to_be_deleted,
                                      const CubitFacetboolOp op);

  CubitStatus dofacetboolean_2bodies(BodySM*& body_in1, 
                             BodySM*& body_in2,
                             BodySM*& body_out,
                             bool keep_old,
                             bool& intersection_found,
                             const CubitFacetboolOp op);

CubitStatus FB_imprint_with_curves(BodySM*& body_in,
                             BodySM*& body_out,                             
                             bool keep_old);
  
  CubitStatus dofacetboolean_2bodies_imprint(BodySM*& body_in1, 
                             BodySM*& body_in2,
                             BodySM*& body_out1,
                             BodySM*& body_out2,                             
                             bool keep_old);

  CubitStatus make_FB_edge_list(DLIList<Curve*> &ref_edge_list);

  void get_edge_list_bbox(CubitBox& edge_list_bbox);

private:

  STD(vector)<FB_Edge*> FB_imprint_edges;
  STD(vector)<FB_Coord*> FB_imprint_edge_coords;
  STD(vector)<FSBoundingBox*> FB_imprint_edge_bboxes;
  
  CubitStatus facetbody_to_facetbool(
                               DLIList<OCCSurface*> &facet_surf_list,
                               STD(vector)<double> &body_verts,
                               STD(vector)<int> &body_connections,
                               STD(vector)<int> *f_c_indices,
                               STD(vector)<OCCSurface *>& fsurfarray,
                               STD(vector)<OCCCurve *>& fcurvearray
                               );                       
  int findcurve(OCCCurve *curve, STD(vector)<OCCCurve *>& fcurvearray);
 
  void make_persistents_webcut(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          BodySM *body_out2,
                                          STD(vector)<OCCSurface *>& fsurfarray,
                                          STD(vector)<OCCCurve *>& fcurvearray,
                                          bool *surfs_in_intersection,
                                          bool *surfs_in_subtraction,
                                          bool *curves_in_intersection,
                                          bool  *curves_in_subtraction
                                         );

  void make_persistents_imprint(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          STD(vector)<OCCSurface *>& fsurfarray,
                                          STD(vector)<OCCCurve *>& fcurvearray
                                         );
 
  void make_persistents_boolean(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          STD(vector)<OCCSurface *>& fsurfarray,
                                          STD(vector)<OCCCurve *>& fcurvearray,
                                          bool *surfs_in_intersection,
                                          bool *surfs_in_subtraction,
                                          bool *curves_in_intersection,
                                          bool  *curves_in_subtraction,
                                          const CubitFacetboolOp op,
                                          bool body_1
                                         );

  void make_persistent_curves(DLIList<OCCCurve*> fcurvelist,
                              STD(vector)<OCCCurve *>& fcurvearray,
                              int n,
                              int which_parent = 0);

  void make_persistent_surfaces(DLIList<OCCSurface*> fsurfaceslist,
                                STD(vector)<OCCSurface *>& fsurfarray,
                                int n,
                                int which_parent = 0);

  int find_coord(double xx, double yy, double zz);

  FSBoundingBox* make_edge_bounding_box(int v0, int v1);
  
    //Separate out the different "lumps" in a given body.  Generally,
    // lumps are disjoint volumes.  They will be converted to separate
    // volumes later in the code.
  CubitStatus separate_lumps( BodySM *body_ptr, bool is_sheet_body );

    //First call separate_lumps to find the disjoing lumps in the body,
    // and then separates the lumps into different bodies.  
  CubitStatus separate_shells_into_bodies( BodySM *body_ptr, 
                                           bool is_sheet_body,
                                           DLIList<BodySM*> &new_bodies );

};

#endif
