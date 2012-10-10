#include "FacetboolInterface.hpp"
#include "FacetDataUtil.hpp"
#include "FBDataUtil.hpp"
#include "FBIntersect.hpp"
#include "IntegerHash.hpp"
#include "FacetBody.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "FacetSurface.hpp"
#include "FBDataUtil.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "FacetModifyEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetCurve.hpp"
#include "CubitFacetEdge.hpp"
#include "TDFacetboolData.hpp"
#include "CurveFacetEvalTool.hpp"
#include "FacetPoint.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "CubitBox.hpp"
#include "FBImprint.hpp"
#include "GfxDebug.hpp"
#include "ChollaEngine.hpp"

FacetboolInterface::FacetboolInterface()
{

}

FacetboolInterface::~FacetboolInterface()
{

}

//===============================================================================
// Function   : webcut_FB
// Member Type: PUBLIC
// Description: webcut a FacetBody with a FacetBool cutter
// Author     : John Fowler
// Date       : 02/04
//===============================================================================
CubitStatus FacetboolInterface::webcut_FB(BodySM *bodysm_ptr,
                                         std::vector<double>& cutter_verts,
                                         std::vector<int>& cutter_connections,
                                         bool cutter_is_plane,
                                         CubitBoolean delete_bodies,
                                         CubitBoolean &intersects,
                                         DLIList<BodySM*>& results_list)
                                         
{
  intersects = CUBIT_TRUE;
  DLIList<FacetSurface*> facet_surf_list;
  std::vector<double> body_verts;
  std::vector<int> body_connections, newbodyfacets;
  std::vector<int> newcutterfacets;
  std::vector<int> f_c_indices1;
  std::vector<FacetSurface *> fsurfarray;
  std::vector<FacetCurve *> fcurvearray;
  CubitStatus status;
  int mydebug = 0;
    FacetModifyEngine *fme = FacetModifyEngine::instance();
    FacetBody *facet_body_ptr = CAST_TO(bodysm_ptr, FacetBody);

    facet_body_ptr->get_surfaces(facet_surf_list);

    status = facetbody_to_facetbool(facet_surf_list,body_verts,
                         body_connections,&f_c_indices1,
                         fsurfarray,fcurvearray);
    if( status == CUBIT_FAILURE )
      return status;

    FBIntersect *intersector = new FBIntersect();
    intersector->set_classify_flag(true);
    if ( cutter_is_plane == true ) intersector->set_body2_planar(); 

    status = intersector->intersect(body_verts,body_connections,
                                    cutter_verts,cutter_connections,
                                    newbodyfacets,newcutterfacets,
                                    &f_c_indices1,
                                    0);
    if( status == CUBIT_FAILURE )
      return status;

    bool *surfs_in_intersection, *surfs_in_subtraction;
    bool *curves_in_intersection, *curves_in_subtraction;
    
    surfs_in_intersection = new bool[1+fsurfarray.size()];
    surfs_in_subtraction = new bool[1+fsurfarray.size()];
    curves_in_intersection = new bool[1+fcurvearray.size()];
    curves_in_subtraction = new bool[1+fcurvearray.size()];
    unsigned int k;
    for ( k = 1; k < 1 + fsurfarray.size(); k++ ) 
        surfs_in_subtraction[k] = surfs_in_intersection[k] = false;
    for ( k = 1; k < 1 + fcurvearray.size(); k++ ) 
        curves_in_subtraction[k] = curves_in_intersection[k] = false;

    status = intersector->get_persistent_entity_info(surfs_in_intersection,
                                 curves_in_intersection,surfs_in_subtraction,
                                 curves_in_subtraction,CUBIT_FB_INTERSECTION,1);
    if( status == CUBIT_FAILURE )
      return status;

    facet_surf_list.reset();

    std::vector<int> surfindex, surfindex2;  
    std::vector<int> curveindex, curveindex2;  
    std::vector<double> vertsout, vertsout2;
    std::vector<int> coordsout, coordsout2;
    
    status = intersector->gather_by_boolean(vertsout,coordsout,
                                     &surfindex,&curveindex,0,
                                     CUBIT_FB_INTERSECTION);  
    if( status == CUBIT_FAILURE )
      return status;
    status = intersector->gather_by_boolean(vertsout2,coordsout2,
                                     &surfindex2,&curveindex2,0,
                                     CUBIT_FB_SUBTRACTION);
    if( status == CUBIT_FAILURE )
      return status;

    //  If there were no intersections
    if ( (vertsout.size()  == 0) || (coordsout.size() == 0) ||
         (vertsout2.size() == 0) || (coordsout2.size() == 0) )
    {
        intersects = CUBIT_FALSE;
        delete intersector;
        delete [] surfs_in_intersection; delete [] surfs_in_subtraction;
        delete [] curves_in_intersection; delete [] curves_in_subtraction;
        return CUBIT_SUCCESS;
    }

    CubitPoint *new_point;
    std::vector<CubitPoint *> points;

    for ( k = 0; k < vertsout.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout[k],
                                                     vertsout[k+1],
                                                     vertsout[k+2] );
      points.push_back(new_point);
    } 

    DLIList <CubitFacet *>facet_list;
    DLIList <CubitPoint *>point_list;
    CubitFacet *facet_ptr;
    if(mydebug)
      GfxDebug::clear();
    for ( k = 0; k < coordsout.size(); k += 3 ) 
    {    
      facet_ptr = new CubitFacetData( points[coordsout[k]],
                                      points[coordsout[k+1]],
                                      points[coordsout[k+2]] ); 
      int cptr[3];
      cptr[0] = curveindex[k];
      cptr[1] = curveindex[k+1];
      cptr[2] = curveindex[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex[k/3],cptr[0],cptr[1],cptr[2],false,
              facet_ptr->is_backwards());
      if(mydebug){
        facet_ptr->debug_draw(CUBIT_YELLOW);
        CubitVector tmp1 = points[coordsout[k]]->coordinates();
        CubitVector tmp2 = points[coordsout[k+1]]->coordinates();
        CubitVector tmp3 = points[coordsout[k+2]]->coordinates();
        
        if(cptr[0]){
          GfxDebug::draw_line(tmp1, tmp2, CUBIT_WHITE);
        }
        if(cptr[1]){
          GfxDebug::draw_line(tmp2, tmp3, CUBIT_GREEN);
        }
        if(cptr[2]){
          GfxDebug::draw_line(tmp3, tmp1, CUBIT_BLUE);
        }
      }
      facet_list.append( facet_ptr );     
    }
    if(mydebug){
      GfxDebug::mouse_xforms();
    }
    points.clear(); //  clear out the points vector since we are through with it.

    double feature_angle;
    int interp_order;
    CubitBoolean smooth_non_manifold, split_surfaces;
    BodySM *body_ptr, *body_ptr2;
      //determine if original body is a sheet body
    Body *tmp_body = CAST_TO(bodysm_ptr->topology_entity(), Body);
    bool is_sheet_body = false;
    if( tmp_body->is_sheet_body() ) 
      is_sheet_body = true;

      //mbrewer:  This is not the best solution.  For now, we are using
      // a no feature angle unless we are working with a sheet-body.  For
      // sheet-bodies, however, we are not marking the vertices that
      // bound curves.  Therefore, if we do not use a feature angle, curves
      // are not split correctly.  The correct fix would be to mark
      // the vertices and maintain them through the webcut.
    if( is_sheet_body ) 
      feature_angle = 135.0;
    else
      feature_angle = 0.0;
    interp_order = 0;
    smooth_non_manifold = CUBIT_TRUE;
    split_surfaces = CUBIT_FALSE;

    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             feature_angle,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         feature_angle, interp_order,
                                         body_ptr);
        if( status == CUBIT_FAILURE )
          return status;

        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }
    DLIList<BodySM*> new_bodies;

    status = separate_shells_into_bodies( body_ptr, is_sheet_body, new_bodies );
    if( status == CUBIT_FAILURE )
      return status;
    for(k=new_bodies.size(); k--;)
      results_list.append( new_bodies.get_and_step() );  
      
    vertsout.clear();
    coordsout.clear();                                 

    facet_list.clean_out();    
    for ( k = 0; k < vertsout2.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout2[k],
                                                   vertsout2[k+1],
                                                   vertsout2[k+2] );
      points.push_back(new_point);
    }  
    for ( k = 0; k < coordsout2.size(); k += 3 ) {    
      facet_ptr = new CubitFacetData( points[coordsout2[k]],
                                      points[coordsout2[k+1]],
                                      points[coordsout2[k+2]] );
      int cptr[3];
      cptr[0] = curveindex2[k];
      cptr[1] = curveindex2[k+1];
      cptr[2] = curveindex2[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex2[k/3],cptr[0],cptr[1],cptr[2],false,
              facet_ptr->is_backwards());
      if(mydebug){
        CubitVector tmp1 = points[coordsout2[k]]->coordinates();
        CubitVector tmp2 = points[coordsout2[k+1]]->coordinates();
        CubitVector tmp3 = points[coordsout2[k+2]]->coordinates();
        if(cptr[0]){
          GfxDebug::draw_line(tmp1, tmp2, CUBIT_WHITE);
        }
        if(cptr[1]){
          GfxDebug::draw_line(tmp2, tmp3, CUBIT_GREEN);
        }
        if(cptr[2]){
          GfxDebug::draw_line(tmp3, tmp1, CUBIT_BLUE);
        }
      }
      facet_list.append( facet_ptr );     
    }
    points.clear(); //  clear out the points vector since we are through with it.

    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             feature_angle,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         feature_angle, interp_order,
                                         body_ptr2);
        if( status == CUBIT_FAILURE )
          return status;
        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }
    new_bodies.clean_out();
    status = separate_shells_into_bodies( body_ptr2, is_sheet_body, new_bodies );
    if( status == CUBIT_FAILURE )
      return status;
    for(k=new_bodies.size(); k--;)
      results_list.append( new_bodies.get_and_step() );  
                               
    if ( delete_bodies == CUBIT_TRUE )
      make_persistents_webcut(bodysm_ptr,body_ptr,body_ptr2,fsurfarray,fcurvearray,
                       surfs_in_intersection,surfs_in_subtraction,
                       curves_in_intersection,curves_in_subtraction);

    vertsout2.clear();
    coordsout2.clear();                                 
                              
    delete intersector;
    delete [] surfs_in_intersection; delete [] surfs_in_subtraction;
    delete [] curves_in_intersection; delete [] curves_in_subtraction;
      
    return CUBIT_SUCCESS;
}

CubitStatus FacetboolInterface::separate_lumps( BodySM *body_ptr, 
                                                bool is_sheet_body)
{
    //get all the shells in 'body_ptr'
    DLIList<FacetShell*> facet_shells;
    FacetBody *facet_body = CAST_TO( body_ptr, FacetBody );
    facet_body->get_shells( facet_shells );

    // some shells here might contain more than one connected
    // patch of surfaces...which is illegal.  Separate them into their own shells.
    bool created_shells = false;
    int k,i;
    for( k=facet_shells.size(); k--; )
    {
      FacetShell *facet_shell = facet_shells.get_and_step();
      DLIList<FacetSurface*> facet_surfs;
      facet_shell->get_surfaces( facet_surfs );
      DLIList<FacetSurface*> connected_patch;
      int max_num_passes = facet_surfs.size();
      int num_passes=0;  //prevents infinite loop
      while( facet_surfs.size() && num_passes < max_num_passes )
      {
        connected_patch.clean_out();
        FacetQueryEngine::instance()->get_connected_patch( facet_surfs, connected_patch );

        if( num_passes == 0 ) 
        {
          int kk;
          for( kk=connected_patch.size(); kk--; )
          {
            FacetSurface *f_surf = connected_patch.get_and_step();
            if( is_sheet_body ) 
              f_surf->set_shell_sense( facet_shell, CUBIT_UNKNOWN );
            else
              f_surf->set_shell_sense( facet_shell, CUBIT_FORWARD );
          }
        }
        else //extract surfaces out of current shell and make them into their own shell
        {

          facet_shell->disconnect_surfaces( connected_patch ); 
         
          DLIList<Surface*> tmp_surfs;
          CAST_LIST( connected_patch, tmp_surfs, Surface ); 
          ShellSM *shellsm_ptr;
          FacetModifyEngine::instance()->make_facet_shell(tmp_surfs, shellsm_ptr);
          if ( shellsm_ptr == NULL )
          {
            PRINT_ERROR("Problems building facet based shell entity.\n");
          }
          else
          {
            FacetShell *tmp_shell = static_cast<FacetShell*>(shellsm_ptr);

            int kk;
            //now for each surface, add sense wrt this new shell
            for( kk=tmp_surfs.size(); kk--;)
            {
              FacetSurface *tmp_facet_surf = CAST_TO(tmp_surfs.get_and_step(), FacetSurface); 
              
              //if it's not a sheet body, tag the surfs wrt the shell, FORWARD
              if( is_sheet_body ) 
                tmp_facet_surf->set_shell_sense( tmp_shell, CUBIT_UNKNOWN );
              else
                tmp_facet_surf->set_shell_sense( tmp_shell, CUBIT_FORWARD );
            }

            created_shells = true;
            FacetLump *facet_lump = static_cast<FacetLump*>( facet_shell->get_lump() ); 
            facet_lump->add_shell( tmp_shell ); 
          }
        }
        num_passes++;
      }
    }
  
    bool created_lumps = false;
    DLIList<FacetShell*> void_shells; 
    int number_regions = 0;
    if( created_shells )
    {
      //determine which shells are regions and which are voids
      //get a point we know is outside the body
      CubitBox bbox;
      CubitVector centroid;
      double vol;
      facet_body->mass_properties( centroid, vol );
      FacetQueryEngine::instance()->create_facet_bounding_box( facet_body, bbox );
      CubitVector point_outside = 2*(bbox.maximum() - centroid);  
      point_outside = centroid + point_outside;
      
      facet_shells.clean_out();
      facet_body->get_shells( facet_shells );
      for( k=facet_shells.size(); k--; )
      {
        FacetShell *facet_shell = facet_shells.get_and_step();
        CubitPointContainment point_cont;
        point_cont = facet_shell->point_containment( point_outside );

        if( point_cont == CUBIT_PNT_OUTSIDE )
          number_regions++;
        else if( point_cont == CUBIT_PNT_INSIDE )
          void_shells.append( facet_shell ); 

        //is a region...if it is the second found region, should be 
        //in its own lump
        if( number_regions>1 && point_cont == CUBIT_PNT_OUTSIDE ) 
        {
          created_lumps = true;
        
          //remove this shell from its current Lump
          FacetLump *tmp_lump = static_cast<FacetLump*>(facet_shell->get_lump());
          tmp_lump->remove_shell( facet_shell ); 

          //make a new lump containing this shell
          Lump *new_lump;
          DLIList<ShellSM*> sm_shells;
          sm_shells.append( static_cast<ShellSM*>(facet_shell) );
          FacetModifyEngine::instance()->make_facet_lump( sm_shells, new_lump );

          //add lump to body
          facet_body->add_lump( static_cast<FacetLump*>(new_lump) );
        }
      }
    }

    if( void_shells.size() && created_lumps )
    {
      //the void shells are still in the original lump.  We need to pair each 
      //up with the region that encloses it 
      
      //make a list of all the region (non-void) shells in the body
      DLIList<FacetShell*> region_shells; 
      facet_body->get_shells( region_shells );
      for(k=void_shells.size(); k--;)
      {
        if( region_shells.move_to( void_shells.get_and_step() ) )
          region_shells.change_to( NULL );
      }
      region_shells.remove_all_with_value( NULL );
    
      //for each void lump...find the region that contains it
      for( k=void_shells.size(); k--; )
      {
        FacetShell *void_shell = void_shells.get_and_step();

        //get a point on the void
        DLIList<FacetSurface*> tmp_surfs;
        void_shell->get_surfaces( tmp_surfs );
        CubitVector point_on_shell;
        
        DLIList<CubitFacet*> facet_list;
        tmp_surfs.get()->tris( facet_list );
        
        point_on_shell = facet_list.get()->center();

        //get the region that contains that point
        for( i=region_shells.size(); i--; )
        {
          FacetShell *region_shell = region_shells.get_and_step();

          tmp_surfs.clean_out();
          region_shell->get_surfaces( tmp_surfs );

          if( region_shell->point_containment( point_on_shell ) 
                                                        == CUBIT_PNT_INSIDE )
          {
            //remove the void shell from it's lump
            FacetLump *tmp_lump = static_cast<FacetLump*>(void_shell->get_lump());
            tmp_lump->remove_shell( void_shell ); 

            //add the void shell to the region's lump
            tmp_lump = static_cast<FacetLump*>(region_shell->get_lump());
            tmp_lump->add_shell( void_shell ); 
            break;
          }
        }
      }
    }
    return CUBIT_SUCCESS;
}

CubitStatus FacetboolInterface::separate_shells_into_bodies(BodySM *body_ptr,
                                                            bool is_sheet_body,
                                                 DLIList<BodySM*> &new_bodies)
{
    //first separate the body into its lumps
  if(!separate_lumps(body_ptr,is_sheet_body))
    return CUBIT_FAILURE;
  
    //split out each lump into it's own body
  DLIList<BodySM*> split_bodies;
  FacetModifyEngine::instance()->split_body( body_ptr, split_bodies);
  new_bodies += split_bodies;
  
  return CUBIT_SUCCESS;
}    

CubitStatus FacetboolInterface::facetbody_to_facetbool(
                               DLIList<FacetSurface*> &facet_surf_list,
                               std::vector<double> &body_verts,
                               std::vector<int> &body_connections,
                               std::vector<int> *f_c_indices,
                               std::vector<FacetSurface *>& fsurfarray,
                               std::vector<FacetCurve *>& fcurvearray
                               )
{
  CubitStatus status;
  int mydebug = 0;
  
  if(mydebug){
    GfxDebug::clear();
  }
  int i, j, k, m, n, vtx[3], hashvalue, *hasharrayptr, hasharraysize, ifoundit;
  IntegerHash *hashobj;
  double xx, yy, zz, xval, yval, zval;
  FacetSurface *facet_surface; 
  DLIList<CubitFacet*> facetlist; 
  DLIList<CubitFacetEdge*> c_edgelist, f_edgelist;
  DLIList<FacetCurve *> curve_list;
  CubitPoint *point;
  CubitFacet *facet;
  FacetCurve *curve;
  CubitFacetEdge *c_edge, *edge;
  int numhashbins, c_index;

  numhashbins = 101;
  status = CUBIT_SUCCESS;
  hashobj = new IntegerHash(numhashbins,20);
                                          
  for ( i = 0; i < facet_surf_list.size(); i++ ) {
    facet_surface = facet_surf_list.get_and_step();
    fsurfarray.push_back(facet_surface); 
    facet_surface->get_my_facetedges(f_edgelist);
    for ( m = f_edgelist.size(); m > 0; m-- ) {
      edge = f_edgelist.get_and_step();
      edge->set_flag(0); //  Initialize edge flags.
      if(mydebug)
        edge->debug_draw(CUBIT_PINK);
    }
    //  Assume that each facet surf occurs only once, so there won't be
    //  duplicates on fsarray.
    facetlist.clean_out();
    facet_surface->tris(facetlist);
    curve_list.clean_out();
    facet_surface->get_curves(curve_list);
    for ( m = 0; m < curve_list.size(); m++ ) {    
      curve = curve_list.get_and_step();
      c_index = findcurve(curve,fcurvearray);
      c_edgelist.clean_out();
      curve->get_facets(c_edgelist);
      for ( n = c_edgelist.size(); n > 0; n-- ) {
        c_edge = c_edgelist.get_and_step();
        c_edge->set_flag(c_index+1);
        if(mydebug)
          c_edge->debug_draw(CUBIT_RED);
        
      }
    }
    int efindex[3];
    for ( j = facetlist.size(); j > 0; j-- ) {
      facet = facetlist.get_and_step();
      for ( k = 0; k < 3; k++ ) {
        point = facet->point(k);
        edge = facet->edge(k);
        if(mydebug){
          if(edge->get_flag()){
            edge->debug_draw(CUBIT_WHITE);
          }
        }
        efindex[k] = edge->get_flag();
        xx = point->x(); 
        yy = point->y();
        zz = point->z();
        //  Get a vertex number.
        hashvalue = FBDataUtil::makeahashvaluefrom_coord(xx,yy,zz,numhashbins);
        hasharrayptr = hashobj->getHashBin(hashvalue,&hasharraysize);
        ifoundit = -1;
        for ( m = 0; m < hasharraysize; m++ ) {
          n = hasharrayptr[m];
          xval = body_verts[3*n];
          yval = body_verts[3*n+1];
          zval = body_verts[3*n+2];
          if ( ( fabs(xval-xx) < 1.e-6 ) && 
          ( fabs(yval-yy) < 1.e-6 ) &&
          ( fabs(zval-zz) < 1.e-6 ) ) {
            ifoundit = n;
            break;
          }
        }
        if ( ifoundit == -1 ) {
          ifoundit = body_verts.size()/3;
          body_verts.push_back(xx);
          body_verts.push_back(yy);
          body_verts.push_back(zz);
          hashobj->addtoHashList(hashvalue,ifoundit);
        }
        vtx[k] = ifoundit;
      }
      body_connections.push_back(vtx[0]);
      body_connections.push_back(vtx[1]);
      body_connections.push_back(vtx[2]);
      f_c_indices->push_back(i+1); //  put the surface index for the facet
      f_c_indices->push_back(efindex[2]); //  edge indices for the facet
      f_c_indices->push_back(efindex[0]);
      f_c_indices->push_back(efindex[1]);            
    }    
  }
  if(mydebug)
    GfxDebug::mouse_xforms();
  delete hashobj;   
  return status;
  
}

//===============================================================================
// Function   : dofacetboolean
// Member Type: PUBLIC
// Description: do facetboolean operations on a list of bodies.
// Author     : John Fowler
// Date       : 02/04
//===============================================================================
CubitStatus FacetboolInterface::dofacetboolean(DLIList<BodySM*>& body_list, 
                             BodySM*& newBody,
                             bool keep_old,
                             const CubitFacetboolOp op)
{
  int k;
  BodySM *body_sm1, *body_sm2, *body_out;
  CubitStatus status = CUBIT_FAILURE;
  bool intersection_found = false;

  body_sm1 = body_list.get_and_step();

  for ( k = body_list.size() - 1; k > 0; k-- ) {
    body_out =0;
    if(!body_sm1)
      body_sm1=body_list.get_and_step();
    else{
      body_sm2 = body_list.get_and_step();  
      status = dofacetboolean_2bodies(body_sm1,body_sm2,body_out,keep_old,
                                      intersection_found,op);
      if ( keep_old == false){
          //if there was an intersection, we want to delete the two
          //original bodies if we are not keeping old.
          //Also, if we we were performing an INTERSECTION we
          //want to delete the two original bodies even if there was
          //no intersection between the two bodies (again if we
          //are not keeping originals).
        if(intersection_found == true || op == CUBIT_FB_INTERSECTION)  { 
          FacetQueryEngine::instance()->
            delete_solid_model_entities(body_sm2);
          FacetQueryEngine::instance()->
            delete_solid_model_entities(body_sm1);
          body_sm2=NULL;
          body_sm1=NULL;
        }
          //if we are not keeping old, we want to delete the second
          //body even if there was no intersection between the two
          //bodies (for a subtration).
        else if(op == CUBIT_FB_SUBTRACTION)  { 
          FacetQueryEngine::instance()->
            delete_solid_model_entities(body_sm2);
          body_sm2=NULL;
        }
      }
        
      if ( body_out ) body_sm1 = body_out;
    }
  }  
  newBody = body_sm1;
  
  return status;  
}

//===============================================================================
// Function   : dofacetboolean_subtract
// Member Type: PUBLIC
// Description: do facetboolean subtract operations on a list of bodies.
// Author     : John Fowler
// Date       : 02/04
//===============================================================================
CubitStatus FacetboolInterface::dofacetboolean_subtract(BodySM*& tool_body, 
                                      DLIList<BodySM*>& from_bodies,
                                      DLIList<BodySM*>& new_bodies,
                                      bool keep_old,
                                      bool* to_be_deleted,
                                      const CubitFacetboolOp op)
{
  int k;
  BodySM *body_sm2, *body_out;
  CubitStatus status = CUBIT_FAILURE;
  bool intersection_found = false;

  for ( k = 0; k < from_bodies.size(); k++ ) {
    body_out = 0;  
    body_sm2 = from_bodies.get_and_step();  
    status = dofacetboolean_2bodies(body_sm2,tool_body,body_out,keep_old,
                            intersection_found,op);
    if ( (keep_old == false) && (intersection_found == true) ) { 
      to_be_deleted[k] = true;
    }

    if ( (status == CUBIT_SUCCESS) && (body_out) ) new_bodies.append(body_out);
  }  
 
  return status;  
}

CubitStatus FacetboolInterface::dofacetboolean_2bodies(BodySM*& body_in1, 
                             BodySM*& body_in2,
                             BodySM*& body_out,
                             bool keep_old,
                             bool& intersection_found,
                             const CubitFacetboolOp op)
{

CubitStatus status;
std::vector<double> body_verts;
std::vector<double> body2_verts;
std::vector<int> body_connections, newbodyfacets;
std::vector<int>  body2_connections, newbody2facets;
std::vector<double> vertsout;
std::vector<int> coordsout;
std::vector<int> f_c_indices1;
std::vector<int> f_c_indices2;
std::vector<FacetSurface *> fsurfarray1, fsurfarray2;
std::vector<FacetCurve *> fcurvearray1, fcurvearray2;
std::vector<int> surfindex;  
std::vector<int> curveindex;
 
 
    DLIList <CubitFacet *>facet_list;
    DLIList <CubitPoint *>point_list;
    CubitPoint *new_point;
    CubitFacet *facet_ptr;
    std::vector<CubitPoint *> points;
    bool *surfs_in_intersection, *surfs_in_subtraction;
    bool *curves_in_intersection, *curves_in_subtraction;
    bool *surfs_in_intersection2, *surfs_in_subtraction2;
    bool *curves_in_intersection2, *curves_in_subtraction2;

  status = CUBIT_SUCCESS;
  intersection_found = true;
  bool is_sheet_body = false;
  Body *tmp_body_1 = CAST_TO(body_in1->topology_entity(), Body);
  Body *tmp_body_2 = CAST_TO(body_in2->topology_entity(), Body);
    //try to figure out if we are using sheet bodies
    //if it is ambiguous, print a warning.
  if(tmp_body_1 && tmp_body_1->is_sheet_body()){
    is_sheet_body = true;
  }
  else if(tmp_body_2 && tmp_body_2->is_sheet_body()){
    is_sheet_body = true;
  }
  if(tmp_body_1 && tmp_body_2 &&
     tmp_body_1->is_sheet_body() != tmp_body_2->is_sheet_body())
  {
    PRINT_WARNING("Geometric boolean requested for a sheet body and a non-sheet body.\n");
  }
  FacetBody *fbody_ptr;
  fbody_ptr = dynamic_cast<FacetBody *>(body_in1);
  DLIList<FacetSurface*> facet_surf_list;
  fbody_ptr->get_surfaces(facet_surf_list);
  
  status = facetbody_to_facetbool(facet_surf_list,body_verts,
                                body_connections,&f_c_indices1,
                                fsurfarray1,fcurvearray1);
  facet_surf_list.clean_out();                                            

    fbody_ptr = dynamic_cast<FacetBody *>(body_in2);
    fbody_ptr->get_surfaces(facet_surf_list);
    status = facetbody_to_facetbool(facet_surf_list,body2_verts,
                                  body2_connections,&f_c_indices2,
                                  fsurfarray2,fcurvearray2);
    facet_surf_list.clean_out();                                    

    FBIntersect *intersector = new FBIntersect();
    intersector->set_classify_flag(true);
    
    status = intersector->intersect(body_verts,body_connections,
                                    body2_verts,body2_connections,
                                    newbodyfacets,newbody2facets,
                                    &f_c_indices1,
                                    &f_c_indices2);

    std::vector<bool> is_body_1;
    status = intersector->gather_by_boolean(vertsout,coordsout,
                                     &surfindex,&curveindex,&is_body_1,op);      

     //  If there were no intersections
    if ( (vertsout.size() == 0) || (coordsout.size() == 0)) {
      delete intersector;
      intersection_found = false;
      return CUBIT_SUCCESS;
    }  
    //  If there was no body_out, we are just checking if there was any
    //  intersection.  If we got this far, there was -- so return.
    if ( &body_out == 0 ) return CUBIT_SUCCESS;
      
    surfs_in_intersection = new bool[1+fsurfarray1.size()];
    surfs_in_subtraction = new bool[1+fsurfarray1.size()];
    curves_in_intersection = new bool[1+fcurvearray1.size()];
    curves_in_subtraction = new bool[1+fcurvearray1.size()];
    surfs_in_intersection2 = new bool[1+fsurfarray2.size()];
    surfs_in_subtraction2 = new bool[1+fsurfarray2.size()];
    curves_in_intersection2 = new bool[1+fcurvearray2.size()];
    curves_in_subtraction2 = new bool[1+fcurvearray2.size()];

    unsigned int k;
    for ( k = 1; k < 1 + fsurfarray1.size(); k++ ) {
        surfs_in_subtraction[k] = surfs_in_intersection[k] = false;
    }
    for ( k = 1; k < 1 + fcurvearray1.size(); k++ ) {
        curves_in_subtraction[k] = curves_in_intersection[k] = false;
    }
    
    for ( k = 1; k < 1 + fsurfarray2.size(); k++ ) {
        surfs_in_subtraction2[k] = surfs_in_intersection2[k] = false;
    }
    for ( k = 1; k < 1 + fcurvearray2.size(); k++ ) {
        curves_in_subtraction2[k] = curves_in_intersection2[k] = false;
    }
     
    status = intersector->get_persistent_entity_info(surfs_in_intersection,
                                 curves_in_intersection,surfs_in_subtraction,
                                 curves_in_subtraction,op,1);
     //  If op == unite, the curves ans surfs that are undamaged are those
     //  for which xxx_in_intereseciton = 0 false and xxx_in_subtraction == true

    status = intersector->get_persistent_entity_info(surfs_in_intersection2,
                                 curves_in_intersection2,surfs_in_subtraction2,
                                 curves_in_subtraction2,op,2);

    delete intersector;

    for ( k = 0; k < vertsout.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout[k],
                                                     vertsout[k+1],
                                                     vertsout[k+2] );
      points.push_back(new_point);
    }  
    for ( k = 0; k < coordsout.size(); k += 3 ) {    
      facet_ptr = new CubitFacetData( points[coordsout[k]],
                                      points[coordsout[k+1]],
                                      points[coordsout[k+2]] );
      int cptr[3];
      cptr[0] = curveindex[k];
      cptr[1] = curveindex[k+1];
      cptr[2] = curveindex[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex[k/3],cptr[0],cptr[1],cptr[2],is_body_1[k/3],
              facet_ptr->is_backwards());                                

      facet_list.append( facet_ptr );     
    }
    points.clear(); //  clear out the points vector since we are through with it.

    FacetModifyEngine *fme = FacetModifyEngine::instance();
    int interp_order;
    CubitBoolean smooth_non_manifold, split_surfaces;
    BodySM *body_ptr;
    interp_order = 0;
    smooth_non_manifold = CUBIT_TRUE;
    split_surfaces = CUBIT_FALSE;
      
    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             -1.0,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         -1.0, interp_order,
                                         body_ptr);
        if( status == CUBIT_FAILURE )
          return status;
        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }

    if ( keep_old == false ) {
      make_persistents_boolean(body_in2,body_ptr,fsurfarray2,fcurvearray2,
                       surfs_in_intersection2,surfs_in_subtraction2,
                       curves_in_intersection2,curves_in_subtraction2,op,false);
      make_persistents_boolean(body_in1,body_ptr,fsurfarray1,fcurvearray1,
                       surfs_in_intersection,surfs_in_subtraction,
                       curves_in_intersection,curves_in_subtraction,op,true);

    }

    body_out = body_ptr;
      //separate the lumps in the "body_out", these will be converted to
      //separate volumes later in the code.
    status=separate_lumps(body_out,is_sheet_body);
    vertsout.clear();
    coordsout.clear();                 
    delete [] surfs_in_intersection; delete [] surfs_in_subtraction;
    delete [] curves_in_intersection; delete [] curves_in_subtraction;
    delete [] surfs_in_intersection2; delete [] surfs_in_subtraction2;
    delete [] curves_in_intersection2; delete [] curves_in_subtraction2;
                                  
  return status; 

}
                           
int FacetboolInterface::findcurve(FacetCurve *curve, 
                                  std::vector<FacetCurve *>& fcurvearray)
{
unsigned int i;

  for ( i = 0; i < fcurvearray.size(); i++ ) {
    if ( fcurvearray[i] == curve ) {
      return i;
    }
  }
  fcurvearray.push_back(curve);
  
  return fcurvearray.size() - 1;
}

void FacetboolInterface::make_persistents_webcut(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          BodySM *body_out2,
                                          std::vector<FacetSurface *>& fsurfarray,
                                          std::vector<FacetCurve *>& fcurvearray,
                                          bool *surfs_in_intersection,
                                          bool *surfs_in_subtraction,
                                          bool *curves_in_intersection,
                                          bool  *curves_in_subtraction
                                         )
{
unsigned int k, n;
DLIList<CubitSimpleAttrib*> csa_list;                              
FacetBody *facet_body_in = CAST_TO(body_in, FacetBody);
FacetBody *facet_body_out1;
FacetBody *facet_body_out2;                               

  facet_body_out1 = CAST_TO(body_out1, FacetBody);
  facet_body_out2 = CAST_TO(body_out2, FacetBody);                               

//  Fix the curves.                      
DLIList<FacetCurve*> fcurvelist, fcurvelist2; 

    facet_body_out1->get_curves(fcurvelist);
    facet_body_out2->get_curves(fcurvelist2);

  for ( n = 1; n < 1 + fcurvearray.size(); n++ ) {
    if ( (curves_in_intersection[n] == true) && (curves_in_subtraction[n] == false) ) {
      make_persistent_curves(fcurvelist,fcurvearray,n);
    } else if ( (curves_in_intersection[n] == false) && 
                (curves_in_subtraction[n] == true) ) {
      make_persistent_curves(fcurvelist2,fcurvearray,n);
    }
  }


//  Fix the surfaces
DLIList<FacetSurface*> fsurfaceslist, fsurfaceslist2;
    facet_body_out1->get_surfaces(fsurfaceslist);
    facet_body_out2->get_surfaces(fsurfaceslist2);

  for ( n = 1; n < 1 + fsurfarray.size(); n++ ) {
    if ( (surfs_in_intersection[n] == true) && (surfs_in_subtraction[n] == false) ) {
      make_persistent_surfaces(fsurfaceslist,fsurfarray,n);                          
    } else if ( (surfs_in_intersection[n] == false) && 
                (surfs_in_subtraction[n] == true) ) {
      make_persistent_surfaces(fsurfaceslist2,fsurfarray,n);                
    }
  }

//  Fix the lumps.
FacetLump *florig, *fl2;
DLIList<FacetLump*> flumplist, flumplist2;
    facet_body_in->get_lumps(flumplist);
      facet_body_out1->get_lumps(flumplist2);
      florig = flumplist.get();
      fl2 = flumplist2.get();
      if ( florig->owner() ) {
        florig->owner()->swap_bridge(florig,fl2,false);
        florig->get_simple_attribute(csa_list);
        for ( k = csa_list.size(); k > 0; k-- ) {
          CubitSimpleAttrib* csa = csa_list.get_and_step();
          fl2->append_simple_attribute_virt(csa);                               
        }
      }
//  Fix the bodies.
    if ( facet_body_in->owner() ) {
      facet_body_in->owner()->swap_bridge(facet_body_in,facet_body_out1,false);
      facet_body_in->get_simple_attribute(csa_list);
      for ( k = csa_list.size(); k > 0; k-- ) {
        CubitSimpleAttrib* csa = csa_list.get_and_step();
        facet_body_out1->append_simple_attribute_virt(csa);                               
      } 
    } 
  
}

void FacetboolInterface::make_persistents_imprint(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          std::vector<FacetSurface *>& fsurfarray,
                                          std::vector<FacetCurve *>& fcurvearray
                                         )
{
DLIList<CubitSimpleAttrib*> csa_list;                              
FacetBody *facet_body_in = CAST_TO(body_in, FacetBody);
FacetBody *facet_body_out1 = CAST_TO(body_out1, FacetBody);                                    
unsigned int n;

//  Fix the curves.                      
    DLIList<FacetCurve*> fcurvelist; 
    facet_body_out1->get_curves(fcurvelist);

  for ( n = 1; n < 1 + fcurvearray.size(); n++ ) {
      make_persistent_curves(fcurvelist,fcurvearray,n,0);
  }


//  Fix the surfaces
    DLIList<FacetSurface*> fsurfaceslist;
    facet_body_out1->get_surfaces(fsurfaceslist);

  for ( n = 1; n < 1 + fsurfarray.size(); n++ ) {
      make_persistent_surfaces(fsurfaceslist,fsurfarray,n,0);                          
  }
 
//  Fix the lumps.
FacetLump *florig, *fl2;
DLIList<FacetLump*> flumplist, flumplist2;
int k;

    facet_body_in->get_lumps(flumplist);
      facet_body_out1->get_lumps(flumplist2);
      florig = flumplist.get();
      fl2 = flumplist2.get();
      if ( florig->owner() ) {
        florig->owner()->swap_bridge(florig,fl2,false);
        florig->get_simple_attribute(csa_list);
        for ( k = csa_list.size(); k > 0; k-- ) {
          CubitSimpleAttrib* csa = csa_list.get_and_step();
          fl2->append_simple_attribute_virt(csa);                               
        }
      }

//  Fix the bodies.
    if ( facet_body_in->owner() ) {
      facet_body_in->owner()->swap_bridge(facet_body_in,facet_body_out1,false);
      facet_body_in->get_simple_attribute(csa_list);
      for ( k = csa_list.size(); k > 0; k-- ) {
        CubitSimpleAttrib* csa = csa_list.get_and_step();
        facet_body_out1->append_simple_attribute_virt(csa);                               
      }
    }    

}

void FacetboolInterface::make_persistents_boolean(BodySM *body_in, 
                                          BodySM *body_out1, 
                                          std::vector<FacetSurface *>& fsurfarray,
                                          std::vector<FacetCurve *>& fcurvearray,
                                          bool *surfs_in_intersection,
                                          bool *surfs_in_subtraction,
                                          bool *curves_in_intersection,
                                          bool  *curves_in_subtraction,
                                          const CubitFacetboolOp op,
                                          bool body_1
                                         )
{
  unsigned int n;
  DLIList<CubitSimpleAttrib*> csa_list;                              
  FacetBody *facet_body_in = CAST_TO(body_in, FacetBody);
  FacetBody *facet_body_out1 = CAST_TO(body_out1, FacetBody);
  int which_parent;
  bool bvalue1 = false, bvalue2 = false;

  if ( body_1 == true ) which_parent = 1;
  else which_parent = 2;

  switch (op) {
  
    case CUBIT_FB_UNION:
      bvalue1 = false;
      bvalue2 = true;
      break;
    case CUBIT_FB_INTERSECTION:
      bvalue2 = false;
      bvalue1 = true;
      break;
    case CUBIT_FB_SUBTRACTION:
      if ( body_1 == true ) {
        bvalue2 = false;
        bvalue1 = true;
      } else {
        bvalue1 = false;
        bvalue2 = true;
      }      
      break;      
  }      
//  Fix the curves.                      
    DLIList<FacetCurve*> fcurvelist; 
    facet_body_out1->get_curves(fcurvelist);

  for ( n = 1; n < 1 + fcurvearray.size(); n++ ) {
    if ( (curves_in_intersection[n] == bvalue1) && 
         (curves_in_subtraction[n] == bvalue2) ) {
      make_persistent_curves(fcurvelist,fcurvearray,n,which_parent);
    }
  }


//  Fix the surfaces
    DLIList<FacetSurface*> fsurfaceslist;
    facet_body_out1->get_surfaces(fsurfaceslist);

  for ( n = 1; n < 1 + fsurfarray.size(); n++ ) {
    if ( (surfs_in_intersection[n] == bvalue1) && 
         (surfs_in_subtraction[n] == bvalue2) ) {
      make_persistent_surfaces(fsurfaceslist,fsurfarray,n,which_parent);                          
    }
  }
 
  if ( body_1 == true ) {
//  Fix the lumps.
FacetLump *florig, *fl2;
DLIList<FacetLump*> flumplist, flumplist2;
int k;

    facet_body_in->get_lumps(flumplist);
      facet_body_out1->get_lumps(flumplist2);
      florig = flumplist.get();
      fl2 = flumplist2.get();
      if ( florig->owner() ) {
        florig->owner()->swap_bridge(florig,fl2,false);
        florig->get_simple_attribute(csa_list);
        for ( k = csa_list.size(); k > 0; k-- ) {
          CubitSimpleAttrib* csa = csa_list.get_and_step();
          fl2->append_simple_attribute_virt(csa);                               
        }
      }

//  Fix the bodies.
    if ( facet_body_in->owner() ) {
      facet_body_in->owner()->swap_bridge(facet_body_in,facet_body_out1,false);
      facet_body_in->get_simple_attribute(csa_list);
      for ( k = csa_list.size(); k > 0; k-- ) {
        CubitSimpleAttrib* csa = csa_list.get_and_step();
        facet_body_out1->append_simple_attribute_virt(csa);                               
      }
    }    
  }
   
}

void FacetboolInterface::make_persistent_curves(DLIList<FacetCurve*> fcurvelist,
                                       std::vector<FacetCurve *>& fcurvearray,
                                       int n,
                                       int which_parent)
{
int k, m, jj;
FacetCurve *fcurveorig, *fcurve2;
int *marked3, index;
bool ifoundit;
FacetPoint *fpointorig, *fpoint2; 
DLIList<FacetPoint*> fpointlist, fpointlist2; 
CubitFacet *cfac;
DLIList<CubitSimpleAttrib*> csa_list;
bool is_from_1, error;
TDFacetboolData* tdf;
 
  fcurveorig = fcurvearray[n-1];
  ifoundit = false;
  for ( k = fcurvelist.size(); k > 0; k-- ) {
    fcurve2 = fcurvelist.get_and_step();
    if ( fcurve2->owner() ) continue;  //  If owner is already set, skip it.
    CurveFacetEvalTool *evaltool = fcurve2->get_eval_tool();
    DLIList<CubitFacetEdge*> edgelist;
    CubitFacetEdge *cfedge;
    evaltool->get_facets(edgelist);
    cfedge = edgelist.get();
    cfac = cfedge->adj_facet(0);
    index = cfac->edge_index(cfedge);
    tdf = TDFacetboolData::get(cfac);
    if ( which_parent ) {
      is_from_1 = tdf->parent_is_body_1();
      if ( ( (is_from_1 == true) && (which_parent == 2) ) ||
           ( (is_from_1 == false) && (which_parent == 1) ) ) continue;
    } 
    marked3 = tdf->get_edge_indices(cfac->is_backwards());
    if ( marked3[(index+1)%3] == n ) ifoundit = true;
    else {
      cfac = cfedge->adj_facet(1);
      if (cfac) {
        index = cfac->edge_index(cfedge);
        tdf = TDFacetboolData::get(cfac);
        if ( which_parent ) {
          is_from_1 = tdf->parent_is_body_1();
          if ( ( (is_from_1 == true) && (which_parent == 2) ) ||
               ( (is_from_1 == false) && (which_parent == 1) ) ) continue;
        }         
        marked3 = tdf->get_edge_indices(cfac->is_backwards());
        if ( marked3[(index+1)%3] == n ) ifoundit = true; 
      }
    }
    if ( ifoundit == true ) {
      fpointlist.clean_out();
      fpointlist2.clean_out();
      fcurveorig->get_points(fpointlist);
      fcurve2->get_points(fpointlist2);
      for ( jj = fpointlist.size(); jj > 0; jj-- ) {
        error = false;
        fpointorig = fpointlist.get_and_step();
        fpoint2 = fpointlist2.get_and_step();
        while ( fpoint2->coordinates() != fpointorig->coordinates() ) {
          if ( fpointlist2.is_at_beginning() == CUBIT_TRUE ) {
//            PRINT_WARNING("Unable to make point on curve persistent.\n");
            error = true;
            break;
          } 
          fpoint2 = fpointlist2.get_and_step();
        }
        if ( error == true ) continue;
        if ( fpointorig->owner() && !(fpoint2->owner()) ) {
          fpointorig->owner()->swap_bridge(fpointorig,fpoint2,false); 
          fpointorig->get_simple_attribute(csa_list);
          for ( m = csa_list.size(); m > 0; m-- ) {
            CubitSimpleAttrib* csa = csa_list.get_and_step();
            fpoint2->append_simple_attribute_virt(csa);                                     
          } 
        }             
      }
      if ( fcurveorig->owner() != 0 ) 
        fcurveorig->owner()->swap_bridge(fcurveorig,fcurve2,false); 
      fcurveorig->get_simple_attribute(csa_list);
      for ( m = csa_list.size(); m > 0; m-- ) {
        CubitSimpleAttrib* csa = csa_list.get_and_step();
        fcurve2->append_simple_attribute_virt(csa);                                     
      }    
      break;
    }
  }
}

void FacetboolInterface::make_persistent_surfaces(DLIList<FacetSurface*> fsurfaceslist,
                            std::vector<FacetSurface *>& fsurfarray,
                            int n,
                            int which_parent)
{
DLIList<CubitSimpleAttrib*> csa_list;                              
int k, m;
FacetSurface *fsorig, *fsurf2;
bool is_from_1;

  fsorig = fsurfarray[n-1];
  for ( k = fsurfaceslist.size(); k > 0; k-- ) {
    fsurf2 = fsurfaceslist.get_and_step();
    DLIList<CubitFacet*> facet_list2;
    fsurf2->tris(facet_list2);
    CubitFacet* facet2 = facet_list2.get();
    TDFacetboolData* tdf = TDFacetboolData::get(facet2);
    if ( which_parent ) {
      is_from_1 = tdf->parent_is_body_1();
      if ( ( (is_from_1 == true) && (which_parent == 2) ) ||
           ( (is_from_1 == false) && (which_parent == 1) ) ) continue;
    }
    int marked2 = tdf->get_surf_index(); 
    if ( marked2 == n ) {
      if ( fsorig->owner() != 0 )
        fsorig->owner()->swap_bridge(fsorig,fsurf2,false); 
      fsorig->get_simple_attribute(csa_list);
      for ( m = csa_list.size(); m > 0; m-- ) {
        CubitSimpleAttrib* csa = csa_list.get_and_step();
        fsurf2->append_simple_attribute_virt(csa);                               
      }    
      break;
    }
  }    
}

CubitStatus FacetboolInterface::FB_imprint_with_curves(BodySM*& body_in,
                             BodySM*& body_out,                             
                             bool keep_old)
{
CubitStatus status;

  FacetModifyEngine *fme = FacetModifyEngine::instance();

  FacetBody *fbody_ptr;
  fbody_ptr = dynamic_cast<FacetBody *>(body_in);
  DLIList<FacetSurface*> facet_surf_list;
  fbody_ptr->get_surfaces(facet_surf_list);
std::vector<double> body_verts;
std::vector<int> body_connections, newbodyfacets;
std::vector<int> f_c_indices;
std::vector<FacetSurface *> fsurfarray;
std::vector<FacetCurve *> fcurvearray;
std::vector<double> vertsout;
std::vector<int> coordsout;
std::vector<int> surfindex;  
std::vector<int> curveindex;  
  
  status = facetbody_to_facetbool(facet_surf_list,body_verts,
                                body_connections,&f_c_indices,
                                fsurfarray,fcurvearray);
  facet_surf_list.clean_out();                                            

  FBImprint *imprinter = new FBImprint();
 
  status = imprinter->imprint_body_curve(body_verts,body_connections,
                                FB_imprint_edge_coords,FB_imprint_edges,
                                FB_imprint_edge_bboxes,&f_c_indices);            

  status = imprinter->update_surfs_and_curves(vertsout,coordsout,
                                          &surfindex,&curveindex);
    DLIList <CubitFacet *>facet_list;
    DLIList <CubitPoint *>point_list;
    CubitPoint *new_point;
    CubitFacet *facet_ptr;
    std::vector<CubitPoint *> points;

unsigned int k;
    for ( k = 0; k < vertsout.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout[k],
                                                     vertsout[k+1],
                                                     vertsout[k+2] );
      points.push_back(new_point);
    }  
    for ( k = 0; k < coordsout.size(); k += 3 ) {        
      facet_ptr = new CubitFacetData( points[coordsout[k]],
                                      points[coordsout[k+1]],
                                      points[coordsout[k+2]] );
      int cptr[3];
      cptr[0] = curveindex[k];
      cptr[1] = curveindex[k+1];
      cptr[2] = curveindex[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex[k/3],cptr[0],cptr[1],cptr[2],false,
              facet_ptr->is_backwards());                                

      facet_list.append( facet_ptr );     
    }
    points.clear(); //  clear out the points vector since we are through with it.
      
    double feature_angle;
    int interp_order;
    CubitBoolean smooth_non_manifold, split_surfaces;
    feature_angle = 135.0;
    interp_order = 0;
    smooth_non_manifold = CUBIT_TRUE;
    split_surfaces = CUBIT_TRUE;
  
    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             feature_angle,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         feature_angle, interp_order,
                                         body_out);
        if( status == CUBIT_FAILURE )
          return status;
        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }

    vertsout.clear();
    coordsout.clear();                                 
    if ( keep_old == false ) {
      make_persistents_imprint(body_in,body_out,fsurfarray,fcurvearray);
    }
    
  delete imprinter;
  
  return status;

}

CubitStatus FacetboolInterface::dofacetboolean_2bodies_imprint(BodySM*& body_in1, 
                             BodySM*& body_in2,
                             BodySM*& body_out1,
                             BodySM*& body_out2,                             
                             bool keep_old)
{
CubitStatus status;
std::vector<double> body_verts;
std::vector<double> body2_verts;
std::vector<int> body_connections, newbodyfacets;
std::vector<int>  body2_connections, newbody2facets;
std::vector<double> vertsout1, vertsout2;
std::vector<int> coordsout1, coordsout2;
std::vector<int> f_c_indices1;
std::vector<int> f_c_indices2;
std::vector<FacetSurface *> fsurfarray1, fsurfarray2;
std::vector<FacetCurve *> fcurvearray1, fcurvearray2;
std::vector<int> surfindex1, surfindex2;  
std::vector<int> curveindex1, curveindex2;  
 
    DLIList <CubitFacet *>facet_list;
    DLIList <CubitPoint *>point_list;
    CubitPoint *new_point;
    CubitFacet *facet_ptr;
    std::vector<CubitPoint *> points;

  FacetModifyEngine *fme = FacetModifyEngine::instance();
  bool is_sheet_body = false;
  Body *tmp_body_1 = CAST_TO(body_in1->topology_entity(), Body);
  Body *tmp_body_2 = CAST_TO(body_in2->topology_entity(), Body);
  if(tmp_body_1 && tmp_body_1->is_sheet_body()){
    is_sheet_body = true;
  }
  else if(tmp_body_2 && tmp_body_2->is_sheet_body()){
    is_sheet_body = true;
  }
  if(tmp_body_1->is_sheet_body() != tmp_body_2->is_sheet_body())
  {
    PRINT_WARNING("Geometric boolean requested for a sheet body and a non-sheet body.\n");
  }
  FacetBody *fbody_ptr;
  fbody_ptr = dynamic_cast<FacetBody *>(body_in1);
  DLIList<FacetSurface*> facet_surf_list;
  fbody_ptr->get_surfaces(facet_surf_list);
  
  status = facetbody_to_facetbool(facet_surf_list,body_verts,
                                body_connections,&f_c_indices1,
                                fsurfarray1,fcurvearray1);
  facet_surf_list.clean_out();                                            

    fbody_ptr = dynamic_cast<FacetBody *>(body_in2);
    fbody_ptr->get_surfaces(facet_surf_list);
    status = facetbody_to_facetbool(facet_surf_list,body2_verts,
                                  body2_connections,&f_c_indices2,
                                  fsurfarray2,fcurvearray2);
    facet_surf_list.clean_out();                                    

    FBIntersect *intersector = new FBIntersect();
    intersector->set_imprint();
    status = intersector->intersect(body_verts,body_connections,
                                    body2_verts,body2_connections,
                                    newbodyfacets,newbody2facets,
                                    &f_c_indices1,
                                    &f_c_indices2);
    status = intersector->update_surfs_and_curves(vertsout1,coordsout1,
                                          &surfindex1,&curveindex1,1);
    status = intersector->update_surfs_and_curves(vertsout2,coordsout2,
                                          &surfindex2,&curveindex2,2);
unsigned int k;
    for ( k = 0; k < vertsout1.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout1[k],
                                                     vertsout1[k+1],
                                                     vertsout1[k+2] );
      points.push_back(new_point);
    }  
    for ( k = 0; k < coordsout1.size(); k += 3 ) {        
      facet_ptr = new CubitFacetData( points[coordsout1[k]],
                                      points[coordsout1[k+1]],
                                      points[coordsout1[k+2]] );
      int cptr[3];
      cptr[0] = curveindex1[k];
      cptr[1] = curveindex1[k+1];
      cptr[2] = curveindex1[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex1[k/3],cptr[0],cptr[1],cptr[2],false,
              facet_ptr->is_backwards());                                

      facet_list.append( facet_ptr );     
    }
    points.clear(); //  clear out the points vector since we are through with it.
      
    double feature_angle;
    int interp_order;
    CubitBoolean smooth_non_manifold, split_surfaces;
    feature_angle = 135.0;
    interp_order = 0;
    smooth_non_manifold = CUBIT_TRUE;
    split_surfaces = CUBIT_FALSE;
  
    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             feature_angle,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         feature_angle, interp_order,
                                         body_out1);
        if( status == CUBIT_FAILURE )
          return status;
        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }
    vertsout1.clear();
    coordsout1.clear();                                 

    facet_list.clean_out();    
    for ( k = 0; k < vertsout2.size(); k += 3 ) {
      new_point = (CubitPoint *) new CubitPointData( vertsout2[k],
                                                   vertsout2[k+1],
                                                   vertsout2[k+2] );
      points.push_back(new_point);
    }  
    for ( k = 0; k < coordsout2.size(); k += 3 ) {    
      facet_ptr = new CubitFacetData( points[coordsout2[k]],
                                      points[coordsout2[k+1]],
                                      points[coordsout2[k+2]] );
      int cptr[3];
      cptr[0] = curveindex2[k];
      cptr[1] = curveindex2[k+1];
      cptr[2] = curveindex2[k+2];
      TDFacetboolData::add_facetbool_facet( facet_ptr );
      TDFacetboolData* td = TDFacetboolData::get(facet_ptr);
      td->set(surfindex2[k/3],cptr[0],cptr[1],cptr[2],false,
              facet_ptr->is_backwards());                              
      facet_list.append( facet_ptr );     
    }
    points.clear(); //  clear out the points vector since we are through with it.
  
    {
        ChollaEngine *cholla_ptr = NULL;

        status = fme->build_cholla_surfaces( facet_list,
                                             point_list,
                                             feature_angle,
                                             interp_order,
                                             smooth_non_manifold,
                                             split_surfaces,
                                             cholla_ptr );
        if( status == CUBIT_FAILURE )
          return status;

        status = fme->finish_facet_Body( cholla_ptr,
                                         NULL,
                                         feature_angle, interp_order,
                                         body_out2);
        if( status == CUBIT_FAILURE )
          return status;
        if ( cholla_ptr )
        {
            cholla_ptr->delete_me();
            delete cholla_ptr;
        }
    }
                               
    if ( keep_old == false ) {
      make_persistents_imprint(body_in1,body_out1,fsurfarray1,fcurvearray1);
      make_persistents_imprint(body_in2,body_out2,fsurfarray2,fcurvearray2);
    }
    
      //separate the lumps in the "body_out1" and "body_out2", these will
      //be converted to separate volumes later in the code.
    status=separate_lumps(body_out1, is_sheet_body);
    status=separate_lumps(body_out2, is_sheet_body);
    vertsout2.clear();
    coordsout2.clear();                                 
                              
    delete intersector;      

  return status;
}

CubitStatus FacetboolInterface::make_FB_edge_list(DLIList<Curve*> &ref_edge_list)
{
CubitStatus success = CUBIT_SUCCESS;
Curve *ref_edge_ptr;
FacetCurve *cfcurve;
CubitPoint  *cfpoint;
DLIList<CubitPoint *> cfpointlist;
int i, j, lastvert, nextvert;
CubitVector coords;

   ref_edge_list.reset();
   for( i = ref_edge_list.size(); i > 0; i-- )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();
      cfcurve = CAST_TO( ref_edge_ptr, FacetCurve);
      cfpointlist.clean_out();
      cfcurve->get_points(cfpointlist);
      cfpoint = cfpointlist.get_and_step();
      coords = cfpoint->coordinates();
      lastvert = find_coord(coords.x(),coords.y(),coords.z());
      for ( j = cfpointlist.size()-1; j > 0; j-- ) {
        cfpoint = cfpointlist.get_and_step();
        coords = cfpoint->coordinates();
        nextvert = find_coord(coords.x(),coords.y(),coords.z());

          FB_Edge *FBedge = new FB_Edge(lastvert,nextvert,0,0,false);
          FSBoundingBox *bb = make_edge_bounding_box(lastvert,nextvert); 
          lastvert = nextvert;         
          FB_imprint_edges.push_back(FBedge);
          FB_imprint_edge_bboxes.push_back(bb);
      }
   }

  return success;  
}

int FacetboolInterface::find_coord(double xx, double yy, double zz)
{
int value;
unsigned int i;

  value = -1;
  for ( i = 0; i < FB_imprint_edge_coords.size(); i++ ) {
    if ( (fabs(FB_imprint_edge_coords[i]->coord[0] - xx) < GEOMETRY_RESABS) &&
         (fabs(FB_imprint_edge_coords[i]->coord[1] - yy) < GEOMETRY_RESABS) &&
         (fabs(FB_imprint_edge_coords[i]->coord[2] - zz) < GEOMETRY_RESABS) ) {
      value = (int)i;
      break;
    }
  }  
  if ( value == -1 ) {
    value = FB_imprint_edge_coords.size();
    FB_Coord *FBcoord = new FB_Coord(xx,yy,zz);
    FB_imprint_edge_coords.push_back(FBcoord);
  }
      
  return value;
}

FSBoundingBox* FacetboolInterface::make_edge_bounding_box(int v0, int v1)
{
FSBoundingBox *bb;
double xmin, ymin, zmin, xmax, ymax, zmax;

  
  xmin = ( FB_imprint_edge_coords[v0]->coord[0] < FB_imprint_edge_coords[v1]->coord[0] ) ?
           FB_imprint_edge_coords[v0]->coord[0] : FB_imprint_edge_coords[v1]->coord[0];
  ymin = ( FB_imprint_edge_coords[v0]->coord[1] < FB_imprint_edge_coords[v1]->coord[1] ) ?
           FB_imprint_edge_coords[v0]->coord[1] : FB_imprint_edge_coords[v1]->coord[1];
  zmin = ( FB_imprint_edge_coords[v0]->coord[2] < FB_imprint_edge_coords[v1]->coord[2] ) ?
           FB_imprint_edge_coords[v0]->coord[2] : FB_imprint_edge_coords[v1]->coord[2];
  xmax = ( FB_imprint_edge_coords[v0]->coord[0] > FB_imprint_edge_coords[v1]->coord[0] ) ?
           FB_imprint_edge_coords[v0]->coord[0] : FB_imprint_edge_coords[v1]->coord[0];
  ymax = ( FB_imprint_edge_coords[v0]->coord[1] > FB_imprint_edge_coords[v1]->coord[1] ) ?
           FB_imprint_edge_coords[v0]->coord[1] : FB_imprint_edge_coords[v1]->coord[1];
  zmax = ( FB_imprint_edge_coords[v0]->coord[2] > FB_imprint_edge_coords[v1]->coord[2] ) ?
           FB_imprint_edge_coords[v0]->coord[2] : FB_imprint_edge_coords[v1]->coord[2];
  if ( (xmax - xmin) < 2.0*GEOMETRY_RESABS ) {
    xmin -= GEOMETRY_RESABS;
    xmax += GEOMETRY_RESABS;
  }
  if ( (ymax - ymin) < 2.0*GEOMETRY_RESABS ) {
    ymin -= GEOMETRY_RESABS;
    ymax += GEOMETRY_RESABS;
  }              
  if ( (zmax - zmin) < 2.0*GEOMETRY_RESABS ) {
    zmin -= GEOMETRY_RESABS;
    zmax += GEOMETRY_RESABS;
  }
  
  bb = new FSBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);

  return bb;

}

void FacetboolInterface::get_edge_list_bbox(CubitBox& edge_list_bbox)
{
unsigned int i;
double min[3],max[3];

  min[0] = min[1] = min[2] = CUBIT_DBL_MAX - 1.0;
  max[0] = max[1] = max[2] = -CUBIT_DBL_MAX;  

  for ( i = 0; i < FB_imprint_edge_coords.size(); i++ ) {
    min[0] = ( min[0] < FB_imprint_edge_coords[i]->coord[0] ) ? min[0] : 
                    FB_imprint_edge_coords[i]->coord[0];
    min[1] = ( min[1] < FB_imprint_edge_coords[i]->coord[1] ) ? min[1] : 
                    FB_imprint_edge_coords[i]->coord[1];
    min[2] = ( min[2] < FB_imprint_edge_coords[i]->coord[2] ) ? min[2] : 
                    FB_imprint_edge_coords[i]->coord[2];                     
    max[0] = ( max[0] > FB_imprint_edge_coords[i]->coord[0] ) ? max[0] : 
                    FB_imprint_edge_coords[i]->coord[0];
    max[1] = ( max[1] > FB_imprint_edge_coords[i]->coord[1] ) ? max[1] : 
                    FB_imprint_edge_coords[i]->coord[1];
    max[2] = ( max[2] > FB_imprint_edge_coords[i]->coord[2] ) ? max[2] : 
                    FB_imprint_edge_coords[i]->coord[2];
  }
  if ( (max[0] - min[0]) < 2.0*GEOMETRY_RESABS ) {
    min[0] -= GEOMETRY_RESABS;
    max[0] += GEOMETRY_RESABS;
  }
  if ( (max[1] - min[1]) < 2.0*GEOMETRY_RESABS ) {
    min[1] -= GEOMETRY_RESABS;
    max[1] += GEOMETRY_RESABS;
  }              
  if ( (max[2] - min[2]) < 2.0*GEOMETRY_RESABS ) {
    min[2] -= GEOMETRY_RESABS;
    max[2] += GEOMETRY_RESABS;
  }   
  CubitBox box(min,max);
  edge_list_bbox = box; 

}



