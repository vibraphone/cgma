//-------------------------------------------------------------------------
// Filename      : FacetSurface.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 06/06/00
//
// Owner         : David R. White
//-------------------------------------------------------------------------
#include "CubitSimpleAttrib.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "FacetSurface.hpp"
#include "RefFace.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetAttrib.hpp"

#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "RefVolume.hpp"
#include "GeometryQueryEngine.hpp"
#include "DLIList.hpp"
#include "FacetEvalTool.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "LoopSM.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "FacetBody.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "FacetLoop.hpp"
#include "FacetCoEdge.hpp"
#include "FacetCurve.hpp"
#include "FacetPoint.hpp"
#include "CubitEvaluator.hpp"
#include "SphereEvaluator.hpp"
#include "CylinderEvaluator.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool. 
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetSurface::FacetSurface(FacetEvalTool *facet_tool,
                           DLIList<ShellSM*> &shellsms,
                           DLIList<LoopSM*> &loopsms)
{
    // Calculate a bounding box if there isn't one already
  facetEvalTool = facet_tool;
    //sense_ = CUBIT_FORWARD;
  myShells += shellsms;
  myLoops += loopsms;
  myShellSense = CUBIT_UNKNOWN;
  myEvaluator = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool. 
//
// Special Notes : Used for save/restore
//
//-------------------------------------------------------------------------
FacetSurface::FacetSurface( const SphereEvaluatorData *sphere_data,
                            FacetEvalTool *facet_tool,
                            DLIList<ShellSM*> &shellsms,
                            DLIList<LoopSM*> &loopsms )
{
  facetEvalTool = facet_tool;
    //sense_ = CUBIT_FORWARD;
  myShells += shellsms;
  myLoops += loopsms;
  myShellSense = CUBIT_UNKNOWN;

  myEvaluator = new SphereEvaluator( sphere_data );
}
//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool. 
//
// Special Notes : Used for save/restore
//
//-------------------------------------------------------------------------
FacetSurface::FacetSurface( const CylinderEvaluatorData *cylinder_data,
                            FacetEvalTool *facet_tool,
                            DLIList<ShellSM*> &shellsms,
                            DLIList<LoopSM*> &loopsms )
{
  facetEvalTool = facet_tool;
    //sense_ = CUBIT_FORWARD;
  myShells += shellsms;
  myLoops += loopsms;
  myShellSense = CUBIT_UNKNOWN;

  myEvaluator = new CylinderEvaluator( cylinder_data );
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool. 
//
// Special Notes : Used for save/restore
//
//-------------------------------------------------------------------------
FacetSurface::FacetSurface(FacetEvalTool *facet_tool,
                           CubitSense sense,
                           CubitSense shell_sense0,
                           CubitBoolean use_facets,
                           DLIList<LoopSM*> &loopsms)
{
    // Calculate a bounding box if there isn't one already
  facetEvalTool = facet_tool;
    //sense_ = CUBIT_FORWARD;
  myLoops += loopsms;
  myShellSense = shell_sense0;
  myEvaluator = NULL;
}
//-------------------------------------------------------------------------
// Purpose       : The default destructor. 
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetSurface::~FacetSurface() 
{
  if ( facetEvalTool )
  {
    delete facetEvalTool;
  }
}


//-------------------------------------------------------------------------
// Purpose       : get the interpolation order of the FacetEvalTool 
//
// Special Notes :
//
//-------------------------------------------------------------------------
int FacetSurface::interp_order() 
{
  assert(facetEvalTool != NULL);
  
  return facetEvalTool->interp_order();
}

//-------------------------------------------------------------------------
// Purpose       : get the min dot of the FacetEvalTool 
//
// Special Notes :
//
//-------------------------------------------------------------------------
double FacetSurface::min_dot() 
{
  assert(facetEvalTool != NULL);
  
  return facetEvalTool->get_min_dot();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void FacetSurface::append_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.append_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void FacetSurface::remove_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.remove_attribute( csa ); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void FacetSurface::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus FacetSurface::get_simple_attribute(DLIList<CubitSimpleAttrib>&
                                                 csa_list)
  { return attribSet.get_attributes(csa_list); }
CubitStatus FacetSurface::get_simple_attribute(const CubitString& name,
                                        DLIList<CubitSimpleAttrib>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus FacetSurface::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }
  
CubitStatus FacetSurface::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }


//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: FacetQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 FacetSurface::get_geometry_query_engine() const
{
   return FacetQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox FacetSurface::bounding_box() const 
{
    if ( myEvaluator )
        return myEvaluator->bounding_box();
    else
        return facetEvalTool->bounding_box();
}


CubitStatus FacetSurface::get_point_normal( CubitVector& location,
                                            CubitVector& normal )
{
  return closest_point( location, NULL, &normal );
}   

CubitStatus FacetSurface::closest_point_uv_guess(  
          CubitVector const& location,
          double& , double& ,
          CubitVector* closest_location,
          CubitVector* unit_normal )
{
  // don't use u and v guesses
 return closest_point(location, closest_location, unit_normal);
}


//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
//-------------------------------------------------------------------------
CubitStatus FacetSurface::closest_point( CubitVector const& location, 
                                         CubitVector* closest_location,
                                         CubitVector* unit_normal_ptr,
                                         CubitVector* curvature1_ptr,
                                         CubitVector* curvature2_ptr)
{
   CubitStatus rv = CUBIT_SUCCESS;
    if ( myEvaluator )
        return myEvaluator->closest_point( location,
                                           closest_location,
                                           unit_normal_ptr,
                                           curvature1_ptr,
                                           curvature2_ptr );

     // Only need to compute the closest location
  if (unit_normal_ptr == NULL &&
      curvature1_ptr  == NULL &&
      curvature2_ptr  == NULL)
  {
    CubitVector temp = location;
    rv=facetEvalTool->closest_point( temp, 
                                  closest_location, NULL );
    if(!rv){
        return rv;
    }
  }
  
    // Need to compute the closest location and the normal, but not the
    // curvatures
  else if ( (unit_normal_ptr != NULL) && 
            (curvature1_ptr == NULL && curvature2_ptr == NULL) )
  {
    CubitVector temp = location;
    rv =facetEvalTool->closest_point( temp, closest_location, unit_normal_ptr);
    if(!rv){
        return rv;
    }
//     if ( get_relative_surface_sense() == CUBIT_REVERSED )
//     {
//       *unit_normal_ptr = -1.0*(*unit_normal_ptr);
//     }
  }
  else
  {
    PRINT_ERROR("Faceted geometry currently doesn't support curvature requests.\n");
    return CUBIT_FAILURE;
  }   
  
  return CUBIT_SUCCESS;
}



CubitStatus FacetSurface::closest_point_along_vector(CubitVector& from_point, 
                                               CubitVector& along_vector,
                                               CubitVector& point_on_surface)
{
  CubitVector other_point = from_point+along_vector;
  DLIList<CubitVector*> intersection_list;
  facetEvalTool->get_intersections( from_point, other_point, intersection_list );
  
  if( intersection_list.size() == 0 )
    return CUBIT_FAILURE;

  if( intersection_list.size() == 1 )
  {
    point_on_surface = *intersection_list.get();
    delete intersection_list.get();
    return CUBIT_SUCCESS;
  }

  //get the closest intersection
  double closest_dist_sq = CUBIT_DBL_MAX;
  for( int k=intersection_list.size(); k--; )
  {
    CubitVector *int_pt = intersection_list.get_and_step();

    double tmp_dist_sq = from_point.distance_between_squared( *int_pt );

    if( tmp_dist_sq < closest_dist_sq )
    {
      point_on_surface = *int_pt;
      closest_dist_sq = tmp_dist_sq;
    }

    delete int_pt;
  }

  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the trimmed surface to the 
//                 input location. 
//
// Special Notes : 
//-------------------------------------------------------------------------
void FacetSurface::closest_point_trimmed( CubitVector from_point, 
                                         CubitVector& point_on_surface)
{
  CubitBoolean on_surf;
  facetEvalTool->closest_point_trimmed( from_point, &point_on_surface, 
                                        on_surf);
  return;
}

//-------------------------------------------------------------------------
// Purpose       : This functions computes the point on the surface that is 
//                 closest to the input location and then calculates the 
//                 magnitudes of the principal curvatures at this (possibly, 
//                 new) point on the surface. Specifying the RefVolume for 
//                 reference is optional.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus FacetSurface::principal_curvatures(
  CubitVector const& location, 
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
  if ( myEvaluator )
  {
      return myEvaluator->principal_curvatures( location,
                                                 curvature_1,
                                                 curvature_2,
                                                 closest_location );
  }

  PRINT_ERROR("Faceted geometry currently does not support curvature requests.\n");
  return CUBIT_FAILURE;
}



CubitStatus FacetSurface::evaluate( double u, double v,
                               CubitVector *position,                                   
                               CubitVector *normal,
                               CubitVector *curvature1,
                               CubitVector *curvature2 )
{
  return CUBIT_FAILURE;
}


//-------------------------------------------------------------------------
// Purpose       : Given values of the two parameters, get the position.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitVector FacetSurface::position_from_u_v (double u, double v)
{
  if ( myEvaluator )
      return myEvaluator->position_from_u_v( u, v );

  PRINT_ERROR("Faceted geometry currently does not support u-v parameterization.\n");
  return CubitVector (0,0,0);
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the {u, v} coordinates of the point 
//                 on the Surface closest to the input point (specified in 
//                 global space). The closest_location is also returned.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus FacetSurface::u_v_from_position( CubitVector const& location,
                                             double& u,
                                             double& v,
                                             CubitVector* closest_location )
{
  if ( myEvaluator )
      return myEvaluator->u_v_from_position( location, u, v, closest_location );

  PRINT_ERROR("Faceted geometry currently does not support u-v parameterization.\n");
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the Facet surface is periodic. Not
//                 available yet.
//                 
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_periodic() 
{
  if ( myEvaluator )
      return myEvaluator->is_periodic();

    //I'm sure we can do something later to calculate this.  But for now...
  //PRINT_ERROR("FacetSurface::is_periodic not implemented yet\n");
  return CUBIT_FALSE;
  
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_periodic_in_U( double& period ) 
{
  if ( myEvaluator )
      return myEvaluator->is_periodic_in_U( period );

  //PRINT_ERROR("FacetSurface::is_periodic_in_U not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_periodic_in_V( double& period ) 
{
  if ( myEvaluator )
      return myEvaluator->is_periodic_in_V( period );

  //PRINT_ERROR("FacetSurface::is_periodic_in_V not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction.  Not available yet.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_singular_in_U( double )
{
  if ( myEvaluator )
      return myEvaluator->is_singular_in_U();

  //PRINT_ERROR("FacetSurface::is_singular_in_U not implemented yet\n");
  return CUBIT_FALSE;
}  

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction.  Not available yet.
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_singular_in_V( double )
{
  if ( myEvaluator )
      return myEvaluator->is_singular_in_V();

  //PRINT_ERROR("FacetSurface::is_singular_in_V not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the U parameter.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_closed_in_U()
{
  if ( myEvaluator )
      return myEvaluator->is_closed_in_U();

  //PRINT_ERROR("FacetSurface::is_closed_in_U not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the V parameter.
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_closed_in_V()
{
  if ( myEvaluator )
      return myEvaluator->is_closed_in_V();

  //PRINT_ERROR("FacetSurface::is_closed_in_V not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Calculates the derivitives at a given parameter location.
//
//-------------------------------------------------------------------------
CubitStatus FacetSurface::uv_derivitives( double ,
                                          double ,
                                          CubitVector &,
                                          CubitVector & )
{
  PRINT_ERROR("Derivitives on a surface not supported on faceted geometry yet.\n");
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the surface is parametrically defined.
//                 Hopefully later this will be available.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_parametric() 
{
  if ( myEvaluator )
  {
      return myEvaluator->is_parametric();
  }

  //PRINT_ERROR("FacetSurface::is_parametric not implemented\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in U, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::get_param_range_U( double& lower_bound,
                                             double& upper_bound )
{
  if ( myEvaluator )
      return myEvaluator->get_param_range_U( lower_bound, upper_bound );

  lower_bound = 0.0;
  upper_bound = 1.0;
  //PRINT_ERROR("FacetSurface::get_param_range_U not implemented\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in V, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::get_param_range_V( double& lower_bound,
                                             double& upper_bound )
{
  if ( myEvaluator )
      return myEvaluator->get_param_range_V( lower_bound, upper_bound );

  lower_bound = 0.0;
  upper_bound = 1.0;
  //PRINT_ERROR("FacetSurface::get_param_range_V not implemented\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns a surface type ID
//
//-------------------------------------------------------------------------
GeometryType FacetSurface::geometry_type()
{
    if ( is_flat() )
    {
        return PLANE_SURFACE_TYPE;
    }
    else if ( is_spherical() )
    {
        return SPHERE_SURFACE_TYPE;
    }
    else if ( is_conical() )
    {
        return CONE_SURFACE_TYPE;
    }
    else
    {
        return FACET_SURFACE_TYPE;
    }
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Surface
//
//-------------------------------------------------------------------------
double FacetSurface::measure() 
{
  return facetEvalTool->area();
}

//-------------------------------------------------------------------------
// Purpose       : Updates the (cached) area of the Surface so that if 
//                 measure() is called, the correct value is returned.
//-------------------------------------------------------------------------
void FacetSurface::update_measurement() 
{
  facetEvalTool->calculate_area();
}

//-------------------------------------------------------------------------
// Purpose       : Returns whether the facet surface is completely flat
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_flat() 
{
  return (facetEvalTool->is_flat() == 1) ? CUBIT_TRUE : CUBIT_FALSE; 
}

//-------------------------------------------------------------------------
// Purpose       : Returns whether the facet surface is spherical
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_spherical()
{
    if ( myEvaluator && myEvaluator->ask_type() == SPHERE_SURFACE_TYPE )
    {
        return CUBIT_TRUE;
    }
    return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns whether the facet surface is spherical
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_conical()
{
    if ( myEvaluator && myEvaluator->ask_type() == CONE_SURFACE_TYPE )
    {
        return CUBIT_TRUE;
    }
    return CUBIT_FALSE;
}

const CubitEvaluatorData *FacetSurface::evaluator_data( void )
{
    if ( myEvaluator )
    {
        return myEvaluator->evaluator_data();
    }
    return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying surface.
//
//-------------------------------------------------------------------------
CubitBoolean FacetSurface::is_position_on( CubitVector &test_position )
{
  CubitVector new_point;
  CubitStatus stat = closest_point(test_position, &new_point, NULL,NULL,NULL);
  if ( !stat )
    return CUBIT_FALSE;
  CubitVector result_vec = test_position - new_point;
  if ( result_vec.length_squared() < GEOMETRY_RESABS )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}
CubitPointContainment FacetSurface::point_containment( const CubitVector &/*point*/ )
{
   return CUBIT_PNT_UNKNOWN;
}
CubitPointContainment FacetSurface::point_containment( double /*u_param*/, 
                                                       double /*v_param*/ )
{
  return CUBIT_PNT_UNKNOWN; 
}
//CubitPointContainment FacetSurface::point_containment( const CubitVector &/*point*/, 
//                                                       double /*u_param*/,
//                                                       double /*v_param*/ )
//{
//   return CUBIT_PNT_UNKNOWN;
//}
CubitSense FacetSurface::get_geometry_sense()
{
    //not sure if this is right for the facet surface...
  CubitSense sense = CUBIT_FORWARD;//get_relative_surface_sense();
  return sense;
}
// CubitSense FacetSurface::get_relative_surface_sense()
// {
//     //not sure if this is right for the facet surface...
//   return sense_;
// }

/*

void FacetSurface::bodysms(DLIList<BodySM*> &bodies) 
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    myShells.get_and_step()->bodysms(bodies);
  }
}

void FacetSurface::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    myShells.get_and_step()->lumps(lumps);
  }
}

void FacetSurface::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    shellsms.append_unique(myShells.get_and_step());
  }
}

void FacetSurface::surfaces(DLIList<Surface*> &surfaces)
{
  surfaces.append_unique(this);
}

void FacetSurface::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    loopsms.append_unique(myLoops.get_and_step());
  }
}
void FacetSurface::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->coedgesms(coedgesms);
  }
}

void FacetSurface::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->curves(curves);
  }
}
void FacetSurface::points(DLIList<Point*> &points)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->points(points);
  }
}
*/


void FacetSurface::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { CAST_LIST_TO_PARENT( myShells, parents ); }
void FacetSurface::get_children_virt( DLIList<TopologyBridge*>& children )
  { CAST_LIST_TO_PARENT( myLoops, children ); }



CubitStatus FacetSurface::get_my_facets(DLIList<CubitFacet*> &facet_list,
                                        DLIList<CubitPoint*>& point_list)
{
  facetEvalTool->get_facets(facet_list);
  facetEvalTool->get_points(point_list);
  return CUBIT_SUCCESS;
}

void FacetSurface::tris(DLIList<CubitFacet *> &facet_list)
{
	facetEvalTool->get_facets(facet_list);
}

void FacetSurface::get_my_points(DLIList<CubitPoint*> &point_list)
{
  facetEvalTool->get_points(point_list);
}

void FacetSurface::get_my_facetedges(DLIList<CubitFacetEdge*> &edge_list)
{
  facetEvalTool->get_edges(edge_list);
}

void FacetSurface::get_shell_sense( CubitSense &sense0)
{
  sense0 = myShellSense;
}


// return the sense with respect to the given shell
CubitSense FacetSurface::get_shell_sense( ShellSM* shell_ptr ) const
{
    // work around non-constness of DLIList functions
  FacetSurface* nonconst = const_cast<FacetSurface*>(this);
  
  int idx = nonconst->myShells.get_index();
  if(idx > 0){
    PRINT_ERROR("Muliple shells attached to surface.\n");
    return CUBIT_UNKNOWN;
  }
  ShellSM *ashell = myShells.get();
  if (ashell == shell_ptr)
  {
    return myShellSense;
  }
//   else
//   {
//     nonconst->myShells.step();
//     idx = nonconst->myShells.get_index();
//     ashell = myShells.get();
//     if (ashell == shell_ptr)
//     {
//       return myShellSense[idx];
//     }
//   }
  return CUBIT_UNKNOWN;
}

// set the sense of the surface with respect to the shell
void FacetSurface::set_shell_sense( FacetShell *facet_shell, 
                                    CubitSense thesense )
{
//    if(thesense == CUBIT_REVERSED){
//        PRINT_INFO("should not do this.");
//    }
  int idx = myShells.get_index();
  if(idx > 0){
    PRINT_ERROR("Multiple shells attached to a single surface.\n");
    return;
  }
  ShellSM *shell_ptr = (ShellSM *)facet_shell;
  ShellSM *ashell = myShells.get();
  if (ashell == shell_ptr)
  {
    myShellSense = thesense;
  }
//   else
//   {
//     myShells.step();
//     idx = myShells.get_index();
//     ashell = myShells.get();
//     if (ashell == shell_ptr)
//     {
//       myShellSense[idx] = thesense;
//     }
//   }
}

//----------------------------------------------------------------
// Function: copy_facets
// Description: copy the points and facets
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetSurface::copy_facets(DLIList<CubitFacet*>&copy_facet_list,
                                      DLIList<CubitPoint*>&copy_point_list)
{
  if (!facetEvalTool)
  {
    PRINT_ERROR("Couldn't copy facets.");
    return CUBIT_FAILURE;
  }
  int ii;
  DLIList<CubitFacet*>facet_list;
  DLIList<CubitPoint*>point_list;
  facetEvalTool->get_facets( facet_list );
  facetEvalTool->get_points( point_list );
  CubitPoint **point_array = new CubitPoint* [point_list.size()];

  //- copy the points

  point_list.reset();
  CubitPoint *new_point, *the_point;
  for(ii=0; ii<point_list.size(); ii++)
  {
    the_point = point_list.get_and_step();
    new_point = new CubitPointData( the_point->coordinates() );
    the_point->marked( ii );
    copy_point_list.append( new_point );
    point_array[ii] = new_point;
  }

  //- copy the facets

  int jj, idx;
  CubitFacet *new_facet, *the_facet;
  CubitPoint *points[3];
  for (ii=0; ii<facet_list.size(); ii++)
  {
    the_facet = facet_list.get_and_step();
    for (jj=0; jj<3; jj++)
    {
      idx = the_facet->point(jj)->marked();
      points[jj] = point_array[idx];
    }
    new_facet = new CubitFacetData( points[0], points[1], points[2] );
    copy_facet_list.append( new_facet );
  }

  delete [] point_array;

  return CUBIT_SUCCESS;
}

void FacetSurface::get_bodies( DLIList<FacetBody*>& result_list )
{
  DLIList<FacetLump*> lump_list;
  get_lumps( lump_list );
  lump_list.reset();
  for ( int i = lump_list.size(); i--; )
  {
    FacetLump* lump = lump_list.get_and_step();
    FacetBody* body = dynamic_cast<FacetBody*>(lump->get_body());
    if (body)
      result_list.append_unique(body);
  }
}

void FacetSurface::get_lumps( DLIList<FacetLump*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    FacetShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    FacetLump* lump = dynamic_cast<FacetLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void FacetSurface::get_shells( DLIList<FacetShell*>& result_list )
{
  myShells.reset();
  for ( int i = 0; i < myShells.size(); i++ )
    if ( FacetShell* shell = dynamic_cast<FacetShell*>(myShells.next(i)) )
      result_list.append(shell);
}


void FacetSurface::get_loops( DLIList<FacetLoop*>& result_list )
{
  myLoops.reset();
  for ( int i = 0; i < myLoops.size(); i++ )
    if ( FacetLoop* loop = dynamic_cast<FacetLoop*>(myLoops.next(i)) )
      result_list.append(loop);
}

void FacetSurface::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  DLIList<FacetLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = 0; i < loop_list.size(); i++ )
    loop_list.next(i)->get_coedges( result_list );
}

void FacetSurface::get_curves( DLIList<FacetCurve*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetCurve* curve = dynamic_cast<FacetCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
}


//-------------------------------------------------------------------------
// Purpose       : Remove Shell from shell list
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus FacetSurface::remove_shell(FacetShell* shell)
{
    // Something strange here -- A DLIList of Shells and a 
    // two-element array for senses?  Try to keep the senses
    // intact anyway...
  myShells.reset();
  if (myShells.get() == shell)
    myShellSense = CUBIT_UNKNOWN;
  
  if (!myShells.move_to(shell))
    return CUBIT_FAILURE;
  
  myShells.remove();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
void FacetSurface::disconnect_all_loops()
{
  myLoops.reset();
  for (int i = myLoops.size(); i--; )
  {
    LoopSM* sm_ptr = myLoops.get_and_step();
    FacetLoop* loop = dynamic_cast<FacetLoop*>(sm_ptr);
    if (loop)
    {
      assert(loop->get_surface() == this);
      loop->remove_surface();
    }
  }
  myLoops.clean_out();
}

void FacetSurface::add_transformation( CubitTransformMatrix &tfmat )
{
    if ( myEvaluator )
        myEvaluator->add_transformation( tfmat );
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
void  FacetSurface::reverse_sense()
{
  facetEvalTool->reverse_facets();
  myLoops.reset();
  int i,j;
  FacetLoop* this_loop;
  LoopSM* this_loop_sm;
  DLIList<FacetCoEdge *> this_coedge_list;
  for(i=0;i<myLoops.size();i++){
    this_loop_sm= myLoops.get_and_step();
    this_loop = dynamic_cast<FacetLoop*>(this_loop_sm);
    if(!this_loop){
      PRINT_ERROR("Unexpected null pointer for loop.\n");
      return;
    }
   this_loop->reverse();
   this_coedge_list.clean_out();
   this_loop->get_coedges(this_coedge_list);
   for(j=0; j<this_coedge_list.size(); j++){
     this_coedge_list.get_and_step()->reverse_sense();
   }
  }
  
  //sense_ = CubitUtil::opposite_sense( sense_ );
  myShellSense = CubitUtil::opposite_sense( myShellSense );
  //myShellSense[1] = CubitUtil::opposite_sense( myShellSense[1] );
}

CubitStatus FacetSurface::get_projected_distance_on_surface( CubitVector *pos1,
                                                             CubitVector *pos2, 
                                                             double &distance )
{
  return CUBIT_FAILURE;
}

CubitStatus FacetSurface::get_sphere_params
(
  CubitVector &center,
  double &radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine sphere parameters for FacetSurfaces.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetSurface::get_cone_params
(
  CubitVector &center,
  CubitVector &normal,
  CubitVector &major_axis,
  double &radius_ratio,
  double &sine_angle,
  double &cos_angle
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine cone parameters for FacetSurfaces.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetSurface::get_torus_params
(
  CubitVector &center,
  CubitVector &normal,
  double &major_radius,
  double &minor_radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine torus parameters for FacetSurfaces.\n");
  return CUBIT_FAILURE;
}


CubitStatus FacetSurface::get_nurb_params
(
  bool &rational,
  int &degree_u,
  int &degree_v,
  int &num_cntrl_pts_u,
  int &num_cntrl_pts_v,
  DLIList<CubitVector> &cntrl_pts,
  DLIList<double> &cntrl_pt_weights,
  DLIList<double> &u_knots,
  DLIList<double> &v_knots
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine nurbs parameters for FacetSurface.\n");
  return CUBIT_FAILURE;
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
