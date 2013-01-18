//------------------------------------------------------------------
//Class: VirtualImprintTool
//Description: Provides functionality and interface to imprinting
//             topology through virtual geometry and faceted booleans.
//Author:  David R. White
//Date: 2/15/2002
//------------------------------------------------------------------

#include "VirtualImprintTool.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "RTree.hpp"
#include "ImprintBoundaryTool.hpp"
#include "AppUtil.hpp"

CubitBoolean VirtualImprintTool::useRealIntersection = CUBIT_FALSE;

//-------------------------------------------------------------------
// Constructor
//-------------------------------------------------------------------
VirtualImprintTool::VirtualImprintTool()
{
}
//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
VirtualImprintTool::~VirtualImprintTool()
{}

//--------------------------------------------------------------------
// virtual_imprint
//  Does the imprinting between two surfaces...
//--------------------------------------------------------------------
CubitStatus VirtualImprintTool::virtual_imprint(RefFace *ref_face1,
                                                RefFace *ref_face2,
                                                DLIList <RefFace*> &results,
                                                double feature_size,
                                                CubitBoolean &curves_modified)
{
  if ( useRealIntersection )
  {
/*    SurfaceImprintTool imp_tool(feature_size);
    CubitStatus stat = imp_tool.imprint_surface(ref_face1, ref_face2, results);
    if (stat != CUBIT_SUCCESS )
      return stat;
    curves_modified = imp_tool.curves_modified();
    return CUBIT_SUCCESS;
*/
    return CUBIT_FAILURE;
  }
  else
  {
    ImprintBoundaryTool ib_tool(ref_face1, ref_face2, feature_size);
    CubitStatus stat=ib_tool.imprint(results, CUBIT_FALSE);
    if (stat != CUBIT_SUCCESS )
      return stat;
    if ( ib_tool.modified_bound_1() || ib_tool.modified_bound_2() )
      curves_modified = CUBIT_TRUE;
    return CUBIT_SUCCESS;
  }
}
//--------------------------------------------------------------------
// virtual_imprint
//  Does the imprinting between multiple surfaces
//--------------------------------------------------------------------
CubitStatus VirtualImprintTool::virtual_imprint(DLIList <RefFace*> &input_faces,
                                                DLIList <RefFace*> &surf_results,
                                                double feature_size,
                                                CubitBoolean &curves_modified)
{
  int ii, jj;
    //Do the n^2 thing, imprint every volume against every other volume...
  DLIList <RefFace*> faces_stack = input_faces, ref_faces_close, tmp_results;
  RefFace *curr_face, *canidate_face;
  while ( faces_stack.size() > 0 )
  {
    curr_face = faces_stack.pop();
    curr_face->marked(1);
    ref_faces_close.clean_out();
    CubitBox curr_box = curr_face->bounding_box();
    for ( jj = input_faces.size(); jj > 0; jj-- )
    {
      canidate_face = input_faces.get_and_step();
      if ( canidate_face == curr_face || canidate_face->marked() )
        continue;
      CubitBox canidate_box = canidate_face->bounding_box();
      if ( curr_box.overlap(feature_size, canidate_box) )
        ref_faces_close.append(canidate_face);
    }
    for ( jj = ref_faces_close.size(); jj > 0; jj-- )
    {
        //check for interrupt.
      if ( AppUtil::instance()->interrupt() )
      {
        break;
      }
      
      tmp_results.clean_out();
      canidate_face = ref_faces_close.get_and_step();
      PRINT_INFO("Imprint Surfaces: %d %d\n", curr_face->id(), canidate_face->id() );
      CubitStatus stat = virtual_imprint(curr_face, canidate_face,
                                         tmp_results, feature_size,
                                         curves_modified);
      if ( stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems imprinting Surface %d and %d\n",
                    curr_face->id(), canidate_face->id());
        return stat;
      }
      if ( tmp_results.size() )
      {
        faces_stack.remove(canidate_face);
        input_faces.remove(curr_face);
        input_faces.remove(canidate_face);
        int ll;
        for ( ll = tmp_results.size(); ll > 0; ll-- )
        {
          canidate_face = tmp_results.get_and_step();
          canidate_face->marked(0);
          faces_stack.append(canidate_face);
          input_faces.append(canidate_face);
        }
        surf_results += tmp_results;
        break;
      }
    }
      //check for interrupt.
    if ( AppUtil::instance()->interrupt() )
    {
        //just break out and still clean up marks and stuff...
      break;
    }
  }
  for ( ii = 0; ii < surf_results.size(); ii++ )
    surf_results.get_and_step()->marked(0);
  DLIList <RefFace*> tmp_list;
  for ( ii = 0; ii < surf_results.size(); ii++ )
  {
    curr_face = surf_results.get_and_step();
    if ( curr_face->marked() )
      continue;
    else
    {
      tmp_list.append(curr_face);
      curr_face->marked(1);
    }
  }
  surf_results.clean_out();
  for ( ii = 0; ii < tmp_list.size(); ii++ )
  {
    curr_face = tmp_list.get_and_step();
    curr_face->marked(0);
    surf_results.append(curr_face);
  }
  return CUBIT_SUCCESS;
}


CubitStatus VirtualImprintTool::virtual_imprint(RefVolume *ref_volume1,
                                                RefVolume *ref_volume2,
                                                DLIList <RefVolume*> &vol_results,
                                                double feature_size,
                                                CubitBoolean &curves_modified)
{
  CubitBoolean ref_1_mod=CUBIT_FALSE, ref_2_mod = CUBIT_FALSE;
  curves_modified = CUBIT_FALSE;
    //First make sure the bounding boxes are close...
  CubitBox ref_vol1_box = ref_volume1->bounding_box();
  CubitBox ref_vol2_box = ref_volume2->bounding_box();
  if ( !ref_vol1_box.overlap(feature_size, ref_vol2_box) )
    return CUBIT_SUCCESS;

    //Now get the list of surfaces for imprinting.
  DLIList <RefFace*> ref_faces1, ref_faces2, ref_faces_close;
  DLIList <RefFace*> faces1_stack, tmp_faces1,tmp_faces2, tmp_output;
  DLIList <RefVolume*> tmp_vols;
  RefFace *curr_ref_face, *canidate_ref_face;
  ref_volume1->ref_faces(ref_faces1);
  ref_volume2->ref_faces(ref_faces2);
  int ii, jj;
  CubitStatus stat;
    //loop over the the ref_faces from the first list.
    //Get all the reffaces that are within tolerance of its
    //bounding box.
  for ( ii = ref_faces1.size(); ii > 0; ii-- )
  {
    curr_ref_face = ref_faces1.get_and_step();
    CubitBox curr_box = curr_ref_face->bounding_box();
    ref_faces_close.clean_out();
	//do this everytime to get the most up-to-date ones...
    ref_faces2.clean_out();
    ref_volume2->ref_faces(ref_faces2);
    for ( jj = ref_faces2.size(); jj > 0; jj-- )
    {
      canidate_ref_face = ref_faces2.get_and_step();
      CubitBox canidate_box = canidate_ref_face->bounding_box();
      if ( curr_box.overlap(feature_size, canidate_box) )
        ref_faces_close.append(canidate_ref_face);
    }
      //Now imprint the faces that are close with the curr_ref_face.
    faces1_stack.clean_out();
    faces1_stack.append( curr_ref_face);
    while( faces1_stack.size() && (curr_ref_face = faces1_stack.pop()) != NULL &&
           ref_faces_close.size() > 0 )
    {
      if ( AppUtil::instance()->interrupt() )
      {
          //just break out and still clean up marks and stuff...
        break;
      }
      while (ref_faces_close.size() && (canidate_ref_face = ref_faces_close.pop()) != NULL )
      {
        tmp_output.clean_out();
        stat = virtual_imprint(curr_ref_face,canidate_ref_face,
                               tmp_output, feature_size, curves_modified );
        if ( stat != CUBIT_SUCCESS )
          return stat;
          //check for interrupt.
        if ( AppUtil::instance()->interrupt() )
        {
            //just break out and still clean up marks and stuff...
          break;
        }
        //Now sort the faces in tmp_output into two lists associated
        //with the owning volumes...
        for ( jj = tmp_output.size(); jj > 0;jj-- )
        {
          RefFace *ref_face = tmp_output.get_and_step();
          tmp_vols.clean_out();
          ref_face->ref_volumes(tmp_vols);
          if ( tmp_vols.move_to(ref_volume1) )
          {
            ref_1_mod = CUBIT_TRUE;
            tmp_faces1.append(ref_face);
          }
          else
          {
            ref_2_mod = CUBIT_TRUE;
            tmp_faces2.append(ref_face);
          }
        }
        if ( tmp_faces1.size() > 0 )
          break;
      }
      faces1_stack += tmp_faces1;
      ref_faces_close += tmp_faces2;
      tmp_faces1.clean_out();
      tmp_faces2.clean_out();
    }
    if ( AppUtil::instance()->interrupt() )
    {
        //just break out and still clean up marks and stuff...
      break;
    }
  }
  if ( ref_1_mod )
    vol_results.append(ref_volume1);
  if ( ref_2_mod )
    vol_results.append(ref_volume2);
    
  return CUBIT_SUCCESS;
}
CubitStatus VirtualImprintTool::virtual_imprint(DLIList <RefVolume*> &input_vols,
                                                DLIList <RefVolume*> &vol_results,
                                                double feature_size,
                                                CubitBoolean &curves_modified)
{
  int ii, jj;
    //Do the n^2 thing, imprint every volume against every other volume...
  DLIList <RefVolume*> marked_vols;
  DLIList <RefVolume*> vols_stack = input_vols, ref_vols_close;
  RefVolume *curr_volume, *canidate_volume;
  while ( vols_stack.size() > 0 )
  {
    if ( AppUtil::instance()->interrupt() )
    {
        //just break out and still clean up marks and stuff...
      break;
    }
    curr_volume = vols_stack.pop();
    curr_volume->marked(1);
    marked_vols.append(curr_volume);
    ref_vols_close.clean_out();
    CubitBox curr_box = curr_volume->bounding_box();
    for ( jj = input_vols.size(); jj > 0; jj-- )
    {
      canidate_volume = input_vols.get_and_step();
      if ( canidate_volume == curr_volume || canidate_volume->marked() )
        continue;
      CubitBox canidate_box = canidate_volume->bounding_box();
      if ( curr_box.overlap(feature_size, canidate_box) )
        ref_vols_close.append(canidate_volume);
    }
    for ( jj = ref_vols_close.size(); jj > 0; jj-- )
    {
      if ( AppUtil::instance()->interrupt() )
      {
          //just break out and still clean up marks and stuff...
        break;
      }
      canidate_volume = ref_vols_close.get_and_step();
      PRINT_INFO("Imprint Volumes: %d %d\n", curr_volume->id(), canidate_volume->id() );
      CubitStatus stat = virtual_imprint(curr_volume, canidate_volume,
                                         vol_results, feature_size,
                                         curves_modified);
      if ( stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems imprinting volume %d and %d\n",
                    curr_volume->id(), canidate_volume->id());
        return stat;
      }
    }
  }
  DLIList <RefVolume*> tmp_list;
  for ( ii = 0; ii < vol_results.size(); ii++ )
  {
    curr_volume = vol_results.get_and_step();
    if ( curr_volume->marked() )
      continue;
    else
    {
      tmp_list.append(curr_volume);
      curr_volume->marked(1);
      marked_vols.append(curr_volume);
    }
  }
  vol_results.clean_out();
  for ( ii = 0; ii < tmp_list.size(); ii++ )
  {
    curr_volume = tmp_list.get_and_step();
    vol_results.append(curr_volume);
  }
  for ( ii = 0; ii < marked_vols.size(); ii++ )
    marked_vols.get_and_step()->marked(0);
  return CUBIT_SUCCESS;
}

    
        
                        
    
  
