#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitVector.hpp"

#include "GeometryModifyEngine.hpp"

CubitStatus GeometryModifyEngine::split_free_curve( Curve *curve, 
                                                    DLIList<CubitVector*> &split_locations,
                                                    DLIList<Curve*> &new_curves )
{
  return CUBIT_FAILURE;
}

Curve* GeometryModifyEngine::make_elliptical_Curve( TBPoint const* point1,
                                                    TBPoint const* point2,
                                                    CubitVector &center_point, 
                                                    double start_angle,
                                                    double end_angle,
                                                    CubitSense sense) const
{
  return NULL;
}

CubitStatus GeometryModifyEngine::create_parallelogram_surface( TBPoint *pt1, 
                                                                TBPoint *pt2,
                                                                TBPoint *pt3,
                                                                BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}
                                                                
CubitStatus GeometryModifyEngine::create_rectangle_surface( double width, 
                                                            double height,
                                                            CubitVector plane,
                                                            BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyEngine::create_circle_surface( TBPoint *pt1, 
                                                         CubitVector center_point,
                                                         TBPoint *pt3,
                                                         BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyEngine::create_circle_surface( TBPoint *pt1, 
                                                         TBPoint *pt3,
                                                         CubitVector center_point,
                                                         BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyEngine::create_circle_surface( double radius,  
                                                         CubitVector plane,
                                                         BodySM *&sheet_body) const
{
  return CUBIT_FAILURE; 
}

CubitStatus GeometryModifyEngine::create_ellipse_surface( TBPoint *pt1, 
                                                          TBPoint *pt3,
                                                          CubitVector center_point,
                                                          BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyEngine::create_ellipse_surface( double major_radius, 
                                                          double minor_radius, 
                                                          CubitVector plane, 
                                                          BodySM *&sheet_body) const
{
  return CUBIT_FAILURE;
}

Curve* GeometryModifyEngine::create_curve_helix( CubitVector &location,
                                                   CubitVector &direction,
                                                   CubitVector &start_point,
                                                   double &thread_distance,
                                                   double &angle,
                                                   bool right_handed) const
{
  return NULL;
}

CubitStatus GeometryModifyEngine::sweep_helical(DLIList<GeometryEntity*>& ref_ent_list,
                                                DLIList<BodySM*>& result_body_list,
                                                CubitVector &location,
                                                CubitVector &direction,
                                                double &thread_distance,
                                                double &angle,
                                                bool right_handed,
                                                bool anchor_entity,
                                                bool keep_old ) const
{
 return CUBIT_FAILURE;
}

