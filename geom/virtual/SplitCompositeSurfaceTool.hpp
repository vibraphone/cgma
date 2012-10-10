#ifndef SPLITCOMPOSITESURFACETOOL_HPP
#define SPLITCOMPOSITESURFACETOOL_HPP

#include "CubitDefines.h"
#include "Surface.hpp"

class Surface;
class RefFace;
class CubitVector;
template<class X> class DLIList;

   
class SplitCompositeSurfaceTool
{
  public:
    static SplitCompositeSurfaceTool *instance();
    static void delete_instance() { if(instance_) delete instance_; };
    CubitStatus split_surface( RefFace *ref_face_ptr,
                             DLIList<CubitVector*> &locations,
                             DLIList<DLIList<CubitVector*>*> &vec_lists,
                             CubitBoolean preview_flg = CUBIT_FALSE,
                             CubitBoolean create_ref_edges_flg = CUBIT_FALSE,
                             CubitBoolean clear_previous_previews = CUBIT_TRUE );

     CubitStatus split_surface( DLIList<RefFace*> &ref_face_list,
                                DLIList<CubitVector*> &locations,
                                DLIList<DLIList<DLIList<CubitVector*>*>*> &list_of_vec_lists,
                                CubitBoolean preview_flg = CUBIT_FALSE,
                                CubitBoolean create_ref_edges_flg =CUBIT_FALSE,
                                CubitBoolean clear_previous_previews = CUBIT_TRUE );

  private:
                                                                                
    void find_faces_for_pos(CubitVector &pos, DLIList<Surface*> surf_list,
                        CubitPointContainment &containment,
                        DLIList<Surface*> &out_list);
    void get_additional_split_points(Surface *surf,
                            DLIList<DLIList<CubitVector*>*> &vec_lists);
    static SplitCompositeSurfaceTool *instance_;
    void find_face_with_non_zero_param_dir(DLIList<Surface*> &surf_list, 
                                        CubitVector &dir,
                                        CubitVector &pos,
                                        Surface *&ret_surf,
                                        double &du, double &dv,
                                        double &step);

};

#endif

