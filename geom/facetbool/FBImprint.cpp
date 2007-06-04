/*
 *
 *
 * Copyright (C) 2004 Sandia Corporation.  Under the terms of Contract DE-AC04-94AL85000
 * with Sandia Corporation, the U.S. Government retains certain rights in this software.
 *
 * This file is part of facetbool--contact via cubit@sandia.gov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *
 */

#include <vector>
#include <math.h>
#include "FBImprint.hpp"
#include "FBPolyhedron.hpp"
#include "FBRetriangulate.hpp"
#include "GeometryDefines.h"
#include "CubitMessage.hpp"
#include "FBDataUtil.hpp"
#include <stdio.h>

FBImprint::FBImprint()
{

}

FBImprint::~FBImprint()
{

}  

CubitStatus FBImprint::imprint_body_curve(const std::vector<double>& Bodycoords,
                           const std::vector<int>& Bodyconnections,
                           const std::vector<FB_Coord*>& FB_imprint_edge_coords,
                           const std::vector<FB_Edge*>& FB_imprint_edges,
                           const std::vector<FSBoundingBox*>& FB_imprint_edge_bboxes,
                           std::vector<int>* indices)
{
CubitStatus status;
bool new_edge_created;

  status = CUBIT_SUCCESS;   

  if ( Bodycoords.size()%3 != 0 ) {
    PRINT_ERROR("Bad coordinates for first part fed to FBImprint.\n");
    return CUBIT_FAILURE;
  }
  if ( Bodyconnections.size()%3 != 0 ) {
    PRINT_ERROR("Bad connection list for first part fed to FBImprint.\n");
    return CUBIT_FAILURE;
  }
 
  f_c_indices = indices;
  poly = new FBPolyhedron;
  status = poly->makepoly(Bodycoords,Bodyconnections,f_c_indices);
  
  status = edges_tri_intersect(FB_imprint_edge_coords,FB_imprint_edges,
                               FB_imprint_edge_bboxes,new_edge_created);

  if ( new_edge_created == true ) {
    std::vector<int> newFacets;
    status = poly->retriangulate(newFacets);
/*    
    FILE *out;
    out = fopen("1Qaz.fct","w");
    int numtris, numverts;
    numverts = poly->verts.size();

    int ii;  

    numtris = 0;
    for ( ii = 0; ii < poly->tris.size(); ii++ ) 
      if ( poly->tris[ii]->dudded == false ) numtris++;
    fprintf(out,"%d %d\n",numverts,numtris);
    for ( ii = 0; ii < poly->verts.size(); ii++ ) 
      fprintf(out,"%d %le %le %le\n",ii+1,poly->verts[ii]->coord[0],
           poly->verts[ii]->coord[1],poly->verts[ii]->coord[2]);
    for ( ii = 0; ii < poly->tris.size(); ii++ ) {
      if ( poly->tris[ii]->dudded == false )
        fprintf(out,"%d %d %d %d\n",ii+1,1+poly->tris[ii]->v0,  
           1+poly->tris[ii]->v1,1+poly->tris[ii]->v2);
    }   
    fclose(out); 
*/        
  }  
  
  return status;

}

CubitStatus FBImprint::edges_tri_intersect(const std::vector<FB_Coord*>& FB_imprint_edge_coords,
                           const std::vector<FB_Edge*>& FB_imprint_edges,
                           const std::vector<FSBoundingBox*>& FB_imprint_edge_bboxes,
                           bool &new_edge_created)
{
  CubitStatus status;
  unsigned int i, j;
  int numboxesfound, *boxlist;
  FSBoundingBox* edgebox;
  double edge_dir[3], edge_0[3], edge_1[3], edge_length = 0.0;
  bool big_angle;

  boxlist = new int[poly->tris.size()];

  status = CUBIT_SUCCESS;
  new_edge_created = false;
  for ( i = 0; i < FB_imprint_edge_bboxes.size(); i++ ) {
     edgebox = FB_imprint_edge_bboxes[i];
     if ( (edgebox->xmax < poly->polyxmin) || 
          (edgebox->xmin > poly->polyxmax) ||
          (edgebox->ymax < poly->polyymin) || 
          (edgebox->ymin > poly->polyymax) ||
          (edgebox->zmax < poly->polyzmin) || 
          (edgebox->zmin > poly->polyzmax) ) continue;
     poly->kdtree->box_kdtree_intersect(*edgebox,&numboxesfound,boxlist);
     if ( numboxesfound > 0 ) {  //  Get a unit vector along the edge.
       edge_0[0] = FB_imprint_edge_coords[FB_imprint_edges[i]->v0]->coord[0];
       edge_1[0] = FB_imprint_edge_coords[FB_imprint_edges[i]->v1]->coord[0];
       edge_dir[0] = edge_1[0] - edge_0[0];
       edge_0[1] = FB_imprint_edge_coords[FB_imprint_edges[i]->v0]->coord[1];
       edge_1[1] = FB_imprint_edge_coords[FB_imprint_edges[i]->v1]->coord[1];
       edge_dir[1] = edge_1[1] - edge_0[1];
       edge_0[2] = FB_imprint_edge_coords[FB_imprint_edges[i]->v0]->coord[2];
       edge_1[2] = FB_imprint_edge_coords[FB_imprint_edges[i]->v1]->coord[2];
       edge_dir[2] = edge_1[2] - edge_0[2];
       edge_length = sqrt(edge_dir[0]*edge_dir[0] + edge_dir[1]*edge_dir[1] + edge_dir[2]*edge_dir[2]);
       if ( edge_length < GEOMETRY_RESABS ) continue;
       edge_dir[0] /= edge_length;
       edge_dir[1] /= edge_length;
       edge_dir[2] /= edge_length;             
     }

     for ( j = 0; j < (unsigned int)numboxesfound; j++ ) { 
        FB_Triangle *tri = poly->tris[boxlist[j]];     
        if ( (edgebox->xmax < tri->boundingbox.xmin) || 
             (edgebox->xmin > tri->boundingbox.xmax) ||
             (edgebox->ymax < tri->boundingbox.ymin) || 
             (edgebox->ymin > tri->boundingbox.ymax) ||
             (edgebox->zmax < tri->boundingbox.zmin) || 
             (edgebox->zmin > tri->boundingbox.zmax) ) continue;
//  Try to get a reasonable value for imprint_res.  It should depend on the size of the 
//  triangle and the length of the edge.
        double tri_size = sqrt( (tri->boundingbox.xmax-tri->boundingbox.xmin)*
                                (tri->boundingbox.xmax-tri->boundingbox.xmin) +
                                (tri->boundingbox.ymax-tri->boundingbox.ymin)*
                                (tri->boundingbox.ymax-tri->boundingbox.ymin) +
                                (tri->boundingbox.zmax-tri->boundingbox.zmin)*
                                (tri->boundingbox.zmax-tri->boundingbox.zmin) );                                
        if ( edge_length < tri_size ) imprint_res = 0.01*edge_length; 
        else imprint_res = 0.01*tri_size;          
//  Flag triangles with planes greater than ~15 degree angle wrt the edge.
        big_angle = false;
        if ( fabs(edge_dir[0]*tri->a +edge_dir[1]*tri->b +edge_dir[2]*tri->c) > 0.25 ) 
          big_angle = true;  
        status = single_edge_tri_intersect(edge_0,edge_1,new_edge_created,tri,big_angle);
     }

  }

  return status;
  
}

//  edge_0 and edge_1 are edge end-points; edge_dir is unit vector from 0 to 1                           
CubitStatus FBImprint::single_edge_tri_intersect(double *edge_0,double *edge_1,
                                                 bool &new_edge_created,
                                                 FB_Triangle *tri,
                                                 bool big_angle) 
{
CubitStatus status;
double distance0, distance1;

  status = CUBIT_SUCCESS;
  distance0 = edge_0[0]*tri->a + edge_0[1]*tri->b + edge_0[2]*tri->c + tri->d;                                                                     
  distance1 = edge_1[0]*tri->a + edge_1[1]*tri->b + edge_1[2]*tri->c + tri->d;                                                                     

  //  If both end-points are farther away from the plane of the triangle than
  //  imprint_res and on the same side, there is no intersection.
  if ( ( (distance0 > imprint_res) && (distance1 > imprint_res) ) ||
       ( (distance0 < -imprint_res) && (distance1 < -imprint_res) ) ) 
    return status;

  //  Check the edge-triangle border closest distances.

  double d0[3], d1[3], s, t, sunclipped, tunclipped;
  double closest_dist, tri_pt[3];
  bool parallel0, parallel1, parallel2;
  int numptsfound = 0;
  double edge_intersection_pt[3][2]; 
  int edge_vert_type[2];
  
  //  Test for closest distance from the edge to each of the triangle edges.
  //  In order for an edge to be generated in the triangle, the intersection
  //  parameter tunclipped has to lie between 0 and 1.  If this condition is met,
  //  the next requirement is that the intersection parameter for the test
  //  edge, sunclipped, has to be between 0 and 1 and the intersection distance
  //  has to be less than imprint_res, or the intersection edge endpoint 
  //  distance (distance0 or distance1) has to be less than imprint_res.
  
  d0[0] = edge_1[0] - edge_0[0];
  d0[1] = edge_1[1] - edge_0[1];
  d0[2] = edge_1[2] - edge_0[2];
  d1[0] = poly->verts[tri->v1]->coord[0] - poly->verts[tri->v0]->coord[0];
  d1[1] = poly->verts[tri->v1]->coord[1] - poly->verts[tri->v0]->coord[1];
  d1[2] = poly->verts[tri->v1]->coord[2] - poly->verts[tri->v0]->coord[2];
  tri_pt[0] = poly->verts[tri->v0]->coord[0];
  tri_pt[1] = poly->verts[tri->v0]->coord[1];
  tri_pt[2] = poly->verts[tri->v0]->coord[2];
  
  closest_dist = FBDataUtil::closest_seg_seg_dist(edge_0,d0,tri_pt,d1,&s,&t,
              &sunclipped,&tunclipped,&parallel0);
              
  if ( (tunclipped >= 0.0) && (tunclipped <= 1.0) ) {
    if ( (sunclipped >= 0.0) && (sunclipped <= 1.0) && 
         (closest_dist <  imprint_res) ) {
      edge_intersection_pt[0][numptsfound] = tri_pt[0] + tunclipped*d1[0];
      edge_intersection_pt[1][numptsfound] = tri_pt[1] + tunclipped*d1[1];
      edge_intersection_pt[2][numptsfound] = tri_pt[2] + tunclipped*d1[2];
      if ( tunclipped == 0.0 ) edge_vert_type[numptsfound] = VERTEX_0;
      else if ( tunclipped == 1.0 ) edge_vert_type[numptsfound] = VERTEX_1;
      else edge_vert_type[numptsfound] = EDGE_0;
      numptsfound++;     
    }
    else if ( sunclipped < 0.0 ) {
      if ( fabs(distance0) < imprint_res ) {
        edge_intersection_pt[0][numptsfound] = edge_0[0] - distance0*tri->a;
        edge_intersection_pt[1][numptsfound] = edge_0[1] - distance0*tri->b;
        edge_intersection_pt[2][numptsfound] = edge_0[2] - distance0*tri->c;
        edge_vert_type[numptsfound] = INTERIOR_VERT;
        numptsfound++;
      } 
    }
    else if ( sunclipped > 1.0 ) {
      if ( fabs(distance1) < imprint_res ) {
        edge_intersection_pt[0][numptsfound] = edge_1[0] - distance1*tri->a;
        edge_intersection_pt[1][numptsfound] = edge_1[1] - distance1*tri->b;
        edge_intersection_pt[2][numptsfound] = edge_1[2] - distance1*tri->c;
        edge_vert_type[numptsfound] = INTERIOR_VERT;
        numptsfound++;      
      }    
    }
  }             
              
  d1[0] = poly->verts[tri->v2]->coord[0] - poly->verts[tri->v1]->coord[0];
  d1[1] = poly->verts[tri->v2]->coord[1] - poly->verts[tri->v1]->coord[1];
  d1[2] = poly->verts[tri->v2]->coord[2] - poly->verts[tri->v1]->coord[2];
  tri_pt[0] = poly->verts[tri->v1]->coord[0];
  tri_pt[1] = poly->verts[tri->v1]->coord[1];
  tri_pt[2] = poly->verts[tri->v1]->coord[2];
  
  closest_dist = FBDataUtil::closest_seg_seg_dist(edge_0,d0,tri_pt,d1,&s,&t,
              &sunclipped,&tunclipped,&parallel1);

  if ( (tunclipped >= 0.0) && (tunclipped <= 1.0) ) {
    if ( (sunclipped >= 0.0) && (sunclipped <= 1.0) && 
         (closest_dist <  imprint_res) ) {
      edge_intersection_pt[0][numptsfound] = tri_pt[0] + tunclipped*d1[0];
      edge_intersection_pt[1][numptsfound] = tri_pt[1] + tunclipped*d1[1];
      edge_intersection_pt[2][numptsfound] = tri_pt[2] + tunclipped*d1[2]; 
      if ( tunclipped == 0.0 ) edge_vert_type[numptsfound] = VERTEX_1;
      else if ( tunclipped == 1.0 ) edge_vert_type[numptsfound] = VERTEX_2;
      else edge_vert_type[numptsfound] = EDGE_1;
      numptsfound++;     
    }
    else if ( sunclipped < 0.0 ) {
      if ( fabs(distance0) < imprint_res ) {
        edge_intersection_pt[0][numptsfound] = edge_0[0] - distance0*tri->a;
        edge_intersection_pt[1][numptsfound] = edge_0[1] - distance0*tri->b;
        edge_intersection_pt[2][numptsfound] = edge_0[2] - distance0*tri->c;
        edge_vert_type[numptsfound] = INTERIOR_VERT;
        numptsfound++;
      } 
    }
    else if ( sunclipped > 1.0 ) {
      if ( fabs(distance1) < imprint_res ) {
        edge_intersection_pt[0][numptsfound] = edge_1[0] - distance1*tri->a;
        edge_intersection_pt[1][numptsfound] = edge_1[1] - distance1*tri->b;
        edge_intersection_pt[2][numptsfound] = edge_1[2] - distance1*tri->c;
        edge_vert_type[numptsfound] = INTERIOR_VERT;
        numptsfound++;      
      }    
    }
  }                          
         
  if ( numptsfound < 2 ) {         
    d1[0] = poly->verts[tri->v0]->coord[0] - poly->verts[tri->v2]->coord[0];
    d1[1] = poly->verts[tri->v0]->coord[1] - poly->verts[tri->v2]->coord[1];
    d1[2] = poly->verts[tri->v0]->coord[2] - poly->verts[tri->v2]->coord[2];
    tri_pt[0] = poly->verts[tri->v2]->coord[0];
    tri_pt[1] = poly->verts[tri->v2]->coord[1];
    tri_pt[2] = poly->verts[tri->v2]->coord[2];

    closest_dist = FBDataUtil::closest_seg_seg_dist(edge_0,d0,tri_pt,d1,&s,&t,
                &sunclipped,&tunclipped,&parallel2);

    if ( (tunclipped >= 0.0) && (tunclipped <= 1.0) ) {
      if ( (sunclipped >= 0.0) && (sunclipped <= 1.0) && 
           (closest_dist <  imprint_res) ) {
        edge_intersection_pt[0][numptsfound] = tri_pt[0] + tunclipped*d1[0];
        edge_intersection_pt[1][numptsfound] = tri_pt[1] + tunclipped*d1[1];
        edge_intersection_pt[2][numptsfound] = tri_pt[2] + tunclipped*d1[2]; 
        if ( tunclipped == 0.0 ) edge_vert_type[numptsfound] = VERTEX_2;
        else if ( tunclipped == 1.0 ) edge_vert_type[numptsfound] = VERTEX_0;
        else edge_vert_type[numptsfound] = EDGE_2;
        numptsfound++;     
      }
      else if ( sunclipped < 0.0 ) {
        if ( fabs(distance0) < imprint_res ) {
          edge_intersection_pt[0][numptsfound] = edge_0[0] - distance0*tri->a;
          edge_intersection_pt[1][numptsfound] = edge_0[1] - distance0*tri->b;
          edge_intersection_pt[2][numptsfound] = edge_0[2] - distance0*tri->c;
          edge_vert_type[numptsfound] = INTERIOR_VERT;
          numptsfound++;
        } 
      }
      else if ( sunclipped > 1.0 ) {
        if ( fabs(distance1) < imprint_res ) {
          edge_intersection_pt[0][numptsfound] = edge_1[0] - distance1*tri->a;
          edge_intersection_pt[1][numptsfound] = edge_1[1] - distance1*tri->b;
          edge_intersection_pt[2][numptsfound] = edge_1[2] - distance1*tri->c;
          edge_vert_type[numptsfound] = INTERIOR_VERT;
          numptsfound++;      
        }    
      }
    }
  }  
  FB_Edge *edge;
  int v10, v11;
  bool exists;

  if ( numptsfound == 2 ) {
    
    tri->dudded = true;
    v10 = poly->addavertex(edge_intersection_pt[0][0],
                           edge_intersection_pt[1][0],
                           edge_intersection_pt[2][0]);
    v11 = poly->addavertex(edge_intersection_pt[0][1],
                           edge_intersection_pt[1][1],
                           edge_intersection_pt[2][1]);
                    
    if ( v10 != v11 ) {
      exists = poly->edge_exists_in_tri(*tri,v10,v11);
      if ( exists == false ) {
        new_edge_created = true;
        edge = new FB_Edge(v10,v11,edge_vert_type[0],edge_vert_type[1],true);
        tri->edge_list.push_back(edge);
        if ( poly->edge_exists(v10,v11) == false )
          poly->intersection_edges.push_back(edge);
      }
    }
  } else if ( numptsfound == 1 ) {
  //  Is it on an edge?
    int edge_type = UNKNOWN;
    int vtype1 = UNKNOWN_VERT, vtype2 = UNKNOWN_VERT;
    int v_other1 = UNKNOWN_VERT, v_other2 = UNKNOWN_VERT;
//     edge_type = UNKNOWN;                      
    if ( edge_vert_type[0] == EDGE_0 ) {
        v_other1 = tri->v0;
        v_other2 = tri->v1;
        edge_type = EDGE_0;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_1;    
    } else if ( edge_vert_type[0] == EDGE_1 ) {
        v_other1 = tri->v1;
        v_other2 = tri->v2;
        edge_type = EDGE_1;
        vtype1 = VERTEX_1; 
        vtype2 = VERTEX_2;
    } else if ( edge_vert_type[0] == EDGE_2 ) {
        v_other1 = tri->v0;
        v_other2 = tri->v2;
        edge_type = EDGE_2;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_2;
    }
    if ( edge_type != UNKNOWN ) {
      tri->dudded = true;
      v10 = poly->addavertex(edge_intersection_pt[0][0],
                             edge_intersection_pt[1][0],
                             edge_intersection_pt[2][0]);
                             
        exists = poly->edge_exists_in_tri(*tri,v_other1,v10);
        if ( exists == false ) {    
          new_edge_created = true;
          edge = new FB_Edge(v_other1,v10,vtype1,edge_type,false);
          tri->edge_list.push_back(edge);
          if ( poly->edge_exists(v_other1,v10) == false )
            poly->intersection_edges.push_back(edge);
        }
        exists = poly->edge_exists_in_tri(*tri,v10,v_other2);
        if ( exists == false ) {    
          new_edge_created = true;
          edge = new FB_Edge(v10,v_other2,edge_type,vtype2,false);
          tri->edge_list.push_back(edge);
          if ( poly->edge_exists(v10,v_other2) == false )
            poly->intersection_edges.push_back(edge);
        } 
     
    
     }
  }

  return status;
}                           
                                                                    
CubitStatus FBImprint::update_surfs_and_curves(std::vector<double>& out_coords,
                                    std::vector<int>& out_connections,
                                    std::vector<int> *out_surf_index,
                                    std::vector<int> *out_curve_index)
{
unsigned int i;

  for ( i = 0; i < poly->verts.size(); i++ ) {
    out_coords.push_back(poly->verts[i]->coord[0]);
    out_coords.push_back(poly->verts[i]->coord[1]);
    out_coords.push_back(poly->verts[i]->coord[2]);
  } 
  for ( i = 0; i < poly->tris.size(); i++ ) {
    if ( poly->tris[i]->dudded == true ) continue;
    out_connections.push_back(poly->tris[i]->v0);
    out_connections.push_back(poly->tris[i]->v1);
    out_connections.push_back(poly->tris[i]->v2);
    out_surf_index->push_back(poly->tris[i]->cubitsurfaceindex);
    out_curve_index->push_back(poly->tris[i]->cubitedge0index);
    out_curve_index->push_back(poly->tris[i]->cubitedge1index);
    out_curve_index->push_back(poly->tris[i]->cubitedge2index);
  }
  return CUBIT_SUCCESS;
}
