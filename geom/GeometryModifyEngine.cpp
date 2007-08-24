#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "BodySM.hpp"
#include "Surface.hpp"

#include "GeometryModifyEngine.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeometryQueryTool.hpp"

CubitStatus GeometryModifyEngine::webcut_with_brick( 
                                     DLIList<BodySM*>& webcut_body_list,
                                     const CubitVector &center,
                                     const CubitVector axes[3],
                                     const CubitVector &extension,
                                     DLIList<BodySM*> &results_list,
                                     bool imprint )
{
   // Create the brick to cut with
   BodySM *cutting_tool_ptr = brick( center, axes, extension );
   if( cutting_tool_ptr == NULL )
      return CUBIT_FAILURE;

   CubitStatus stat =
       this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_cylinder(
                                        DLIList<BodySM*> &webcut_body_list,
                                        double radius,
                                        const CubitVector &axis,
                                        const CubitVector &center,
                                        DLIList<BodySM*>& results_list,
                                        bool imprint ) 
{
  
  double max_size =  0.;
  //lets find the distance to the center for each body and take
  //the max.
  double curr;
  CubitVector cent_bod;
  CubitBox bounding_box; 
  BodySM *body_ptr;
  bounding_box = webcut_body_list[0]->bounding_box();
  cent_bod =  bounding_box.center();
  cent_bod = cent_bod - center;
  curr = cent_bod.length();
  if ( curr > max_size )
     max_size = curr;


  for ( int ii = webcut_body_list.size()-1; ii > 0; ii-- )
    {
      body_ptr = webcut_body_list[ii];
      bounding_box |= body_ptr->bounding_box();
      cent_bod = body_ptr->bounding_box().center();
      cent_bod = cent_bod - center;
      curr = cent_bod.length();
      if ( curr > max_size )
	max_size = curr;
    }
  
  curr = bounding_box.diagonal().length();

  if ( curr > max_size )
     max_size = curr;

  double height = 0.0;   
  if ( center.x() > max_size )
    {
      height = 500.0 * center.x();
    }
  else if ( center.y() > max_size )
    {
      height = 500.0 * center.y();
    }
  else if ( center.z() > max_size )
    {
      height = 500.0 * center.z();
    }
  else
    {
      height = 500.0 * max_size;
    }

  //lets make certain we have a valid height..
  if ( height < GEOMETRY_RESABS )
    {
      height = 500.0;
    }

  BodySM *cutting_tool_ptr = cylinder( height, radius, radius, radius );

  if( cutting_tool_ptr == NULL )
    return CUBIT_FAILURE;

  //transform the cyclinder to cernter and axis
  // The current frustum is centered on the z axis.
  CubitVector axis2(0., 0., 1.0 );
  //now find the normal to the current axis and axis we want to be
  //at. This normal is where we will rotate about.
  CubitVector normal_axis = axis2 * axis;
  if ( normal_axis.length() > CUBIT_RESABS )
    {
       //angle in degrees.
       double angle = normal_axis.vector_angle( axis2, axis );
       get_gqe()->rotate(cutting_tool_ptr, normal_axis, angle);
    }
  get_gqe()->translate(cutting_tool_ptr, center);
  
  CubitStatus stat =
    this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

  // Delete the BodySM that was created to be used as a tool
  get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

  return stat;
}
                           
CubitStatus GeometryModifyEngine::webcut_with_planar_sheet(
                                     DLIList<BodySM*>& webcut_body_list,
                                     const CubitVector &center,
                                     const CubitVector axes[2],
				     double width,
				     double height,
                                     DLIList<BodySM*> &results_list,
                                     bool imprint )
{
   // Create the planar sheet to cut with
   CubitVector p1, p2, p3, p4;
 
   // Get the corners of the sheet
   center.next_point( axes[0], width/2.0, p1 );
   p1.next_point( axes[1], -height/2.0, p1 );
   p1.next_point( axes[1], height, p2 );
   p2.next_point( axes[0], -width, p3 );
   p3.next_point( axes[1], -height, p4 );

   BodySM *cutting_tool_ptr = planar_sheet(p1,p2,p3,p4);  
   if( cutting_tool_ptr == NULL )
      return CUBIT_FAILURE;

   CubitStatus stat = 
       this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_sweep_curves_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& point,
                              const CubitVector& sweep_axis,
                              double angle,
                              Surface *stop_surf,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  if(curves.size() == 0 )
    return CUBIT_FAILURE;
 
  DLIList<GeometryEntity*> ref_ent_list;
  for(int i = 0; i < curves.size(); i++)    
    ref_ent_list.append((GeometryEntity*)(curves.get_and_step()));			  
  
  DLIList<BodySM*> swept_bodies;
  CubitStatus stat = sweep_rotational(ref_ent_list,swept_bodies,point,
                                      sweep_axis,angle,0, 0.0, 0,false,false,
                                      false,stop_surf);
  if(stat == CUBIT_FAILURE  && swept_bodies.size() == 0)
    return stat;

  //stitch faces together
  BodySM* cutting_tool_ptr = NULL;
  stat = stitch_surfs(swept_bodies, cutting_tool_ptr);

  if(stat == CUBIT_FAILURE)
    {
      //delete all swept faces
      for(int i = 0; i < swept_bodies.size(); i++)
	get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
      return stat;
    }

  stat = 
       this->webcut(blank_bodies, cutting_tool_ptr, results_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   //delete all swept faces
   for(int i = 0; i < swept_bodies.size(); i++)
      get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
   return stat;
  
}

CubitStatus GeometryModifyEngine::webcut_with_sweep_curves(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& sweep_vector,
                              bool through_all,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  if(curves.size() == 0 )
    return CUBIT_FAILURE;
  
  CubitVector tmp_sweep_vector = sweep_vector;

  //get model bbox info...will scale sweep vector by its diagonal
  //so that we go far enough
  if( through_all || stop_surf )
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.normalize();
      tmp_sweep_vector*=(2*bounding_box.diagonal().length());
    }

  DLIList<GeometryEntity*> ref_ent_list;
  for(int i = 0; i < curves.size(); i++)    
    ref_ent_list.append((GeometryEntity*)(curves.get_and_step()));
  DLIList<BodySM*> swept_bodies;
  CubitStatus stat;

  //see if we're sweeping along a specified curve
  if( curve_to_sweep_along )
    {
      DLIList<Curve*> curves_to_sweep_along;
      curves_to_sweep_along.append(curve_to_sweep_along);
      stat = sweep_along_curve(ref_ent_list, swept_bodies,
			       curves_to_sweep_along,0.0,0,false,stop_surf);

      if (!stat && swept_bodies.size() == 0)
	return stat;
    }

  else
    {
      stat = sweep_translational(ref_ent_list, swept_bodies,
			     tmp_sweep_vector,0.0,0, false, false, stop_surf);

      if (!stat && swept_bodies.size() == 0)
	return stat;
    }

  //stitch faces together
  BodySM* cutting_tool_ptr = NULL;
  stat = stitch_surfs(swept_bodies, cutting_tool_ptr);

  if(stat == CUBIT_FAILURE)
    {
      //delete all swept faces
      for(int i = 0; i < swept_bodies.size(); i++)
	get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
      return stat;
    }

  stat = 
       this->webcut(blank_bodies, cutting_tool_ptr, results_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   //delete all swept faces
   for(int i = 0; i < swept_bodies.size(); i++)
      get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
   return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_curve_loop(
			      DLIList<BodySM*> &webcut_body_list,
                              DLIList<Curve*> &curves,
                              DLIList<BodySM*>& results_list,
                              bool imprint )
{
  if(curves.size() == 0 )
    return CUBIT_FAILURE;

  GeometryType surface_type = PLANE_SURFACE_TYPE;
  CubitStatus stat;

  //make copies of the curves.
  DLIList<Curve*> copied_curves;
  Curve* temp_curve = NULL;
  for (int i = curves.size(); i--;)
    {
      temp_curve = make_Curve(curves.get_and_step());
      if(temp_curve != NULL)
        copied_curves.append(temp_curve);
    }
  
  //make a face out of the curves
  Surface * surf = make_Surface(surface_type, copied_curves, NULL, false );
  if (surf == NULL)
    {
      PRINT_ERROR("webcut tool surface is not created from acis.\n");
      return CUBIT_FAILURE;
    }

  
  //get cutting tool BodySM.
  BodySM* cutting_tool_ptr = make_BodySM(surf);
  assert(cutting_tool_ptr );
    
  stat = 
       this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

   // Delete the BodySM that was created to be used as a tool
   get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_extended_surf(
			      DLIList<BodySM*> &webcut_body_list,
                              Surface *extend_from,
                              DLIList<BodySM*> &results_list,
                              int &num_cut,
                              bool imprint )
{
  assert(extend_from);
  //make the extended face 
  Surface * surf = make_Surface(extend_from, true);
  if (surf == NULL)
    {
      PRINT_ERROR("webcut tool surface is not created from acis.\n");
      return CUBIT_FAILURE;
    }

  //get cutting tool BodySM.
  BodySM* cutting_tool_ptr = surf->bodysm();
  assert(cutting_tool_ptr );
    
  CubitStatus stat;
  stat = 
    this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

  num_cut = results_list.size();

  // Delete the BodySM that was created to be used as a tool
  get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
  return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_sweep_surfaces_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& point,
                              const CubitVector& sweep_vector,
                              double angle,
                              Surface *stop_surf,
                              bool up_to_next,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  if(surfaces.size() == 0 )
    return CUBIT_FAILURE;
 
  DLIList<GeometryEntity*> ref_ent_list;
  Surface * temp_face = NULL;
  for(int i = 0; i < surfaces.size(); i++) 
    {
      //copy the faces before sweep
      temp_face = make_Surface(surfaces.get_and_step());
      if (temp_face)
	ref_ent_list.append((GeometryEntity*)temp_face);		       
    }

  BodySM* to_body = NULL;
  CubitStatus stat = CUBIT_SUCCESS;
  if(up_to_next && blank_bodies.size() > 1) //unite all bland_bodies
    {
       DLIList<BodySM*> newBodies;

       DLIList<BodySM*> copied_bodies;
       for(int i = 0; i < blank_bodies.size(); i++)
         copied_bodies.append(copy_body(blank_bodies.get_and_step()));

       stat = unite(copied_bodies, newBodies);

       if(stat == CUBIT_FAILURE)
	 {
	   PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
	   PRINT_INFO("Try the 'stop surface <id>' option instead\n");
	   return stat;
	 }
       to_body = newBodies.get();
    }

  else if(up_to_next && blank_bodies.size() == 1)
    to_body = copy_body(blank_bodies.get());

  DLIList<BodySM*> swept_bodies;
  stat = sweep_rotational(ref_ent_list,swept_bodies,point,
                                      sweep_vector,angle,0, 0.0,0,false,false,
                                      false,stop_surf, to_body);
  if(stat == CUBIT_FAILURE && swept_bodies.size() == 0)
    {
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }

       return stat;
    }

  //if there are more than 1, unite them all
  DLIList<BodySM*>  newBodies;
  if (swept_bodies.size() > 1)
     stat = unite(swept_bodies, newBodies);
  else
     newBodies = swept_bodies;

  if(stat == CUBIT_FAILURE || newBodies.size()!= 1)
    { 
       PRINT_ERROR("webcut tool body is not created from acis.\n");
       //delete the swept_bodies
       BodySM* tmp_body = NULL;
       for (int i = swept_bodies.size(); i--;)
         {
	   tmp_body= swept_bodies.get_and_step();
           if (tmp_body)
	     get_gqe()->delete_solid_model_entities(tmp_body);
	 }
 
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              get_gqe()->delete_solid_model_entities( (Surface*)temp_entity); 
         }
       return CUBIT_FAILURE;
    }

  BodySM *cutting_tool_ptr = newBodies.get(); 

  stat = 
    this->webcut(blank_bodies, cutting_tool_ptr, results_list, imprint ) ;

  // Delete the BodySM that was created to be used as a tool
  get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

  return stat;
}

CubitStatus GeometryModifyEngine::webcut_with_sweep_surfaces(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& sweep_vector,
                              bool sweep_perp,
                              bool through_all,
                              bool outward,
                              bool up_to_next,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  if(surfaces.size() == 0 )
    return CUBIT_FAILURE;
  
  CubitVector tmp_sweep_vector = sweep_vector;

  //get model bbox info...will scale sweep vector by its diagonal
  //so that we go far enough
  if( through_all || stop_surf || up_to_next)
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.normalize();
      tmp_sweep_vector*=(2*bounding_box.diagonal().length());
    }

  if( sweep_perp == true )
  {
    if( through_all || stop_surf )
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.set(1,0,0);
      tmp_sweep_vector = 2*(bounding_box.diagonal());
    }
  }

  BodySM* to_body = NULL;
  CubitStatus stat = CUBIT_SUCCESS;
  if(up_to_next && blank_bodies.size() > 1) //unite all bland_bodies
    {
       DLIList<BodySM*> newBodies;
       DLIList<BodySM*> copied_bodies;
       for(int i = 0; i < blank_bodies.size(); i++)
	 copied_bodies.append(copy_body(blank_bodies.get_and_step()));

       stat = unite(copied_bodies, newBodies);
       if(stat == CUBIT_FAILURE)
	 {
	   PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
	   PRINT_INFO("Try the 'stop surface <id>' option instead\n");
	   return stat;
	 }
       to_body = newBodies.get();
    }

  else if(up_to_next && blank_bodies.size() == 1)
    to_body = copy_body(blank_bodies.get());

  DLIList<GeometryEntity*> ref_ent_list;
  Surface * temp_face = NULL;
  for(int i = 0; i < surfaces.size(); i++) 
    {
      //copy the faces before sweep
      temp_face = make_Surface(surfaces.get_and_step());
      if (temp_face)
	ref_ent_list.append((GeometryEntity*)temp_face);		       
    }
  
  //Sweep surfaces
  DLIList<BodySM*> swept_bodies;

  //see if we're sweeping along a specified curve
  if( curve_to_sweep_along )
    {
      DLIList<Curve*> curves_to_sweep_along;
      curves_to_sweep_along.append(curve_to_sweep_along);
      stat = sweep_along_curve(ref_ent_list, swept_bodies,
			       curves_to_sweep_along, 0.0,0,false,stop_surf,
			       to_body);   
    }

  else if (sweep_perp )
    stat = sweep_perpendicular(ref_ent_list, swept_bodies,
			       tmp_sweep_vector.length(),0.0,0,outward,false,
			       stop_surf, to_body);
  else    
    stat = sweep_translational(ref_ent_list, swept_bodies,
			       tmp_sweep_vector,0.0,0, false, false, stop_surf,
			       to_body);


  if (stat == CUBIT_FAILURE && swept_bodies.size() == 0)
    {
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
       {
          temp_entity = ref_ent_list.get_and_step();
          if (temp_entity)
            get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
       }
       return stat;
    }

  //if there are more than 1, unite them all
  DLIList<BodySM*>  newBodies;
  if (swept_bodies.size() > 1)
     stat = unite(swept_bodies, newBodies);
  else
     newBodies = swept_bodies;

  if(stat == CUBIT_FAILURE || newBodies.size()!= 1)
    { 
       PRINT_ERROR("the webcut tool body is not created from acis.\n");
       //delete the swept_bodies
       BodySM* tmp_body = NULL;
       for (int i = swept_bodies.size(); i--;)
         {
	   tmp_body= swept_bodies.get_and_step();
           if (tmp_body)
	     get_gqe()->delete_solid_model_entities(tmp_body);
	 }
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }

       return CUBIT_FAILURE;
    }

  BodySM *cutting_tool_ptr = newBodies.get();  

  stat = 
    this->webcut(blank_bodies, cutting_tool_ptr, results_list, imprint ) ;

  // Delete the BodySM that was created to be used as a tool
  get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

  return stat;
}
