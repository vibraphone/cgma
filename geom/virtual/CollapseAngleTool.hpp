//---------------------------------------------------------------------------
//- Filename:       CollapseAngleTool
//- Purpose:  To collapse small angle for preparing for mesh
//-
//- Creator:       Jiangtao Hu
//- Creation date: 02/09/2005
//---------------------------------------------------------------------------- 

#ifndef COLLAPSEANGLETOOL_HPP
#define COLLAPSEANGLETOOL_HPP

#include "DLIList.hpp"
class RefVertex;
class RefEdge;
class RefFace;
class CubitVector;                                                 
template<class X> class DLIList;
   
class CollapseAngleTool
{
public:
  static CollapseAngleTool *instance();
    // return a static instance pointer

  CubitStatus auto_collapse (double length,
                             double angle,
                             CubitBoolean preview,
                             CubitBoolean if_comp_vertex,
                             double max_angle);

  void collapse_one_angle(RefVertex *vex_ptr,
                          RefEdge   *edge_to_remove,
                          RefEdge   *the_other_edge,
                          double     length1,
                          double     length2,
                          CubitBoolean get_position,
                          CubitVector  &position,
                          CubitBoolean preview,
                          CubitBoolean if_comp_vertex,
                          double       angle);

private:
  CollapseAngleTool();
    // Constructor
                                                                                
  ~CollapseAngleTool();
    // Destructor
                                                                                
  static CollapseAngleTool *instance_;
    // the static instance pointer

  DLIList <RefEdge *> partly_drawn_curve;
  CubitStatus partition_curve(RefEdge *&edge,
                              RefVertex *root_vertex,
                              CubitVector  &position,
                              RefVertex *&new_vertex);
                                                                               
  CubitBoolean if_partition_surf(RefVertex *vex_ptr,
                                 RefEdge  *edge_to_remove,
                                 RefFace  *common_face);

  CubitStatus partition_curve(RefEdge *&edge,
                              RefVertex *root_vertex,
                              double  arc_length,
                              RefVertex *&new_vertex);
                                                                                
  CubitStatus partition_surface(RefFace *common_face,
                                RefEdge *& edge,
                                RefVertex * root,
                                RefVertex *vertex1,
                                RefVertex *veretx2,
                                RefFace  *&result_face);
                                                                                
  CubitStatus collapse_angle(RefVertex *vex_ptr,
                             RefEdge   *edge_to_remove,
                             RefEdge   *the_other_edge,
                             double     length1,
                             double     length2,
                             CubitBoolean get_position,
                             CubitVector &position);

  CubitStatus draw_preview(RefVertex *vex_ptr,
                           RefEdge   *edge_to_remove,
                           RefEdge   *the_other_edge,
                           double     length1,
                           double     length2,
                           CubitBoolean get_position,
                           CubitVector &position);
                                                                                
  void  exchange_edges(RefEdge *&edge,
                       RefEdge *&the_other_edge);
                                                                               
  double the_surface_angle(RefEdge *edge,
                           RefVertex *root_vertex,
                           double    length);
                                                                                
  CubitStatus position_from_length(RefEdge *edge,
                                   RefVertex *root_vertex,
                                   double  arc_length,
                                   CubitVector& v_new);
};

#endif

