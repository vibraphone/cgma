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

#ifdef WIN32
#  pragma warning(disable : 4786)
#endif

#include <math.h>
#include <deque>

#include "FBPolyhedron.hpp"
#include "FBStructs.hpp"
#include "FBRetriangulate.hpp"
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include "CubitVector.hpp"
#include "GfxDebug.hpp"
const double BOX_CRACK = 1.e-4;

FBPolyhedron::FBPolyhedron()
{

  polyxmin = polyymin = polyzmin = CUBIT_DBL_MAX;
  polyxmax = polyymax = polyzmax = -polyxmin;
  hashobj = new IntegerHash(NUMHASHBINS,20);
  original_numtris = 0;
  
}

FBPolyhedron::~FBPolyhedron()
{
unsigned int i;

  delete hashobj; 
  delete kdtree;
  for ( i = 0; i < verts.size(); i++ ) {
    delete verts[i];
  } 
  for ( i = 0; i < tris.size(); i++ ) {
    delete tris[i];  
  } 
 
}

CubitStatus FBPolyhedron::makepoly(const std::vector<double>& coords,
                                   const std::vector<int>& connections,
                                   std::vector<int> *f_c_indices)
{
  int hashvalue, parent, cubitfacetindex;
  int cubitedge0index, cubitedge1index, cubitedge2index;
  unsigned int i;
  FB_Coord *mycoord; 
  FB_Triangle *mytri;
  CubitStatus status;
  FSBOXVECTOR boxvector;
  std::vector<int>::iterator dpi;
  status = CUBIT_SUCCESS;
  
   for ( i = 0; i < coords.size(); i += 3 ) {

      mycoord = new FB_Coord(coords[i],coords[i+1],coords[i+2]); 
      hashvalue = makeahashvaluefrom_coord(coords[i],coords[i+1],coords[i+2]);
      hashobj->addtoHashList(hashvalue,verts.size());
      verts.push_back(mycoord);  
   }
   numpts = verts.size();
   parent = -1;
   if ( f_c_indices ){
     dpi = f_c_indices->begin();
   }
   for ( i = 0; i < connections.size(); i += 3 ) {
       //as a safety check make, not only do we see if we have the
       // f_c_indices vector, but make sure we don't go past the end of
       // that vector.
     if ( f_c_indices  && (i<f_c_indices->size())){
       cubitfacetindex = *dpi;
       cubitedge0index = *(dpi+1);
       cubitedge1index = *(dpi+2);
       cubitedge2index = *(dpi+3);
       dpi += 4;    
     }
     else {
       cubitfacetindex = 0;
       cubitedge0index = cubitedge1index = cubitedge2index = 0;
     }
     mytri = new FB_Triangle(connections[i],connections[i+1],connections[i+2],
                             parent,cubitfacetindex,cubitedge0index,
                             cubitedge1index,cubitedge2index);
     make_tri_boundingbox(mytri);

     boxvector.push_back(&mytri->boundingbox);  // for the KdTree

     make_tri_plane_coeffs(mytri);     
    
     tris.push_back(mytri);
     
   }
 
   numtris = original_numtris = tris.size();
   kdtree = new KDTree();
   kdtree->makeKDTree(numtris,boxvector);
   
   return status;        
}
      
bool FBPolyhedron::boxintersection(FBPolyhedron *otherpoly)
{
double xmin, ymin, zmin, xmax, ymax, zmax;

  xmin = otherpoly->polyxmin;
  ymin = otherpoly->polyymin;
  zmin = otherpoly->polyzmin;
  xmax = otherpoly->polyxmax;
  ymax = otherpoly->polyymax;
  zmax = otherpoly->polyzmax;
  
  if ( (polyxmin > xmax) || (polyxmax < xmin) ) return false;
  if ( (polyymin > ymax) || (polyymax < ymin) ) return false;
  if ( (polyzmin > zmax) || (polyzmax < zmin) ) return false;
   
  return true;
}

int FBPolyhedron::makeahashvaluefrom_coord(double x, double y, double z)
{
double mantissasum;

      if ( fabs(x) < 1.e-3 ) x = 0.0;
      if ( fabs(y) < 1.e-3 ) y = 0.0;
      if ( fabs(z) < 1.e-3 ) z = 0.0;
      mantissasum = (int)(10000.0*fabs(x) + 0.5) + 
                    (int)(10000.0*fabs(y) + 0.5) + 
          (int)(10000.0*fabs(z) + 0.5);
      
      return (int)(mantissasum) % NUMHASHBINS;
}

int FBPolyhedron::addavertex(double x, double y, double z)
{
int hashvalue, i, ifoundit;
int *hasharrayptr, hasharraysize;
double xval, yval, zval;
FB_Coord *mycoord, *newcoord;

  hashvalue = makeahashvaluefrom_coord(x,y,z);
  
  hasharrayptr = hashobj->getHashBin(hashvalue,&hasharraysize);
  
  ifoundit = -1;
  for ( i = 0; i < hasharraysize; i++ ) {
    mycoord = verts[hasharrayptr[i]];
    xval = mycoord->coord[0];
    yval = mycoord->coord[1];
    zval = mycoord->coord[2];
    if ( ( fabs(xval-x) < EPSILON ) && 
    ( fabs(yval-y) < EPSILON ) &&
    ( fabs(zval-z) < EPSILON ) ) {
      ifoundit = hasharrayptr[i];
      break;
    }
  }
  if ( ifoundit == -1 ) {
    newcoord = new FB_Coord(x,y,z);
    ifoundit = verts.size();    
    verts.push_back(newcoord);
    hashobj->addtoHashList(hashvalue,ifoundit);
  }
  verts[ifoundit]->is_on_boundary = true;
  return ifoundit;
  
}

void FBPolyhedron::putnewpoints(std::vector<double>& newpoints)
{
double coordinate;
unsigned int i;

//  numpts is the original number of points.  Any new points were added  
//  onto the end of the verts vector.

  if ( verts.size() > numpts ) {
    for ( i = numpts; i < verts.size(); i++ ) {
      coordinate = verts[i]->coord[0];
      newpoints.push_back(coordinate);
      coordinate = verts[i]->coord[1];
      newpoints.push_back(coordinate);
      coordinate = verts[i]->coord[2];
      newpoints.push_back(coordinate);
    }  
  }
}

void FBPolyhedron::putedges(std::vector<int>& newedges)
{
std::deque<unsigned int> edgedeque;
std::multimap<unsigned int,unsigned int>::iterator p3;
unsigned int startedge, startvert, edgnum;
std::multimap<unsigned int,unsigned int>::iterator pub;
unsigned int first, second;
unsigned int i, numedgesadded;
int v0, v1;
bool foundone;

  for ( i = 0; i < intersection_edges.size(); i++ ) {
    v0 = intersection_edges[i]->v0;
    v1 = intersection_edges[i]->v1;

    edgnum = i;
    edgmmap.insert(std::pair<const unsigned int, unsigned int>(v0,edgnum));        
    edgmmap.insert(std::pair<const unsigned int, unsigned int>(v1,edgnum));    
  }
    numedgesadded = 0;
    startedge = 0;

  while ( numedgesadded < intersection_edges.size() ) {
    if ( intersection_edges[startedge]->mark == true ) {
      startedge++;
      continue;
    }
    startvert = v0 = intersection_edges[startedge]->v0;
    v1 = intersection_edges[startedge]->v1;
    intersection_edges[startedge]->mark = true;
    numedgesadded += 1;
    edgedeque.push_back(v0);
    edgedeque.push_back(v1);
 
    foundone = true;

    while ( foundone == true ) {       
      p3 = edgmmap.find(v1);
      pub = edgmmap.upper_bound(v1);

      if ( p3 != edgmmap.end() ) {
        do {
          foundone = false;
     first = p3->first;
     second = p3->second;
     if ( (second == startedge) && (first == startvert) ) {
       break;
     }
     if ( (intersection_edges[second]->v0 != v0) &&
          (intersection_edges[second]->v1 != v0) ) {
       if ( intersection_edges[second]->mark == true ) {
         p3++;
         continue;
       }
       if (intersection_edges[second]->v0 == v1 ) {
         v1 = intersection_edges[second]->v1;
         v0 = intersection_edges[second]->v0;
       } else {
         v1 = intersection_edges[second]->v0;
         v0 = intersection_edges[second]->v1;
       }
            intersection_edges[second]->mark = true;
            numedgesadded += 1;
            edgedeque.push_back(v1);
            foundone = true;
     }
     p3++;
        } while ( (p3 != pub) && (foundone == false) );
        if ( foundone == false ) break;
      } 
    } //                     end while ( foundone == true )
    v1 = startvert;
    foundone = true;
    while ( foundone == true ) {       
        p3 = edgmmap.find(v1);
        pub = edgmmap.upper_bound(v1);

        if ( p3 != edgmmap.end() ) {
          do {
            foundone = false;
       first = p3->first;
       second = p3->second;
       if ( (intersection_edges[second]->v0 == v1) ||
            (intersection_edges[second]->v1 == v1) ) {
         if ( intersection_edges[second]->mark == true ) {
           p3++;
           continue;
         }
         if (intersection_edges[second]->v0 == v1 ) {
           v1 = intersection_edges[second]->v1;
           v0 = intersection_edges[second]->v0;
         } else {
           v1 = intersection_edges[second]->v0;
           v0 = intersection_edges[second]->v1;
         }
              intersection_edges[second]->mark = true;
         numedgesadded += 1;
              edgedeque.push_front(v1);
         foundone = true;
       }
       p3++;
          } while ( (p3 != pub) && (foundone == false) );
          if ( foundone == false ) break;
        } 
    }
  
    unsigned int size = edgedeque.size();

    if ( size > 0 ) {
      newedges.push_back(size);
      for ( i = 0; i < size; i++ ) {
        newedges.push_back(edgedeque[i]);
      }
    }
    edgedeque.clear();
    startedge++;
  } //                     end while ( numedgesadded < intersection_edges.size() )
 
  newedges.push_back(0);


}

bool FBPolyhedron::edge_exists(int v0, int v1) 
{
bool exists = false;
unsigned int i;

  for ( i = 0; i < intersection_edges.size(); i++ ) {
    if ( ( (intersection_edges[i]->v0 == v0) && 
           (intersection_edges[i]->v1 == v1) ) ||
         ( (intersection_edges[i]->v0 == v1) && 
           (intersection_edges[i]->v1 == v0) ) ) {
      exists = true;
      break;
    }   
  }

  return exists;
}

CubitStatus FBPolyhedron::retriangulate(std::vector<int>& newfacets, 
                                        std::vector<int>& newfacetsindex)
{
CubitStatus status;
FBRetriangulate *retriangulater;
unsigned int i;

  status = CUBIT_SUCCESS;
  for ( i = 0; i < tris.size(); i++ ) {
    if ( tris[i]->dudded == true ) {
      tris[i]->parent = (int)i;
      retriangulater = new FBRetriangulate(verts, tris, newfacets, newfacetsindex); 
      status = retriangulater->retriangulate_this_tri(i);      
      delete retriangulater;
    }
  
  }
  
  return status;
}

CubitStatus FBPolyhedron::retriangulate(std::vector<int>& newfacets)
{
CubitStatus status;
FBRetriangulate *retriangulater;
unsigned int i;

  status = CUBIT_SUCCESS;
  for ( i = 0; i < tris.size(); i++ ) {
    if ( tris[i]->dudded == true ) {
      tris[i]->parent = (int)i;
      retriangulater = new FBRetriangulate(verts, tris, newfacets); 
      status = retriangulater->retriangulate_this_tri(i);
      delete retriangulater;
    }
  
  }
  
  return status;
}

bool FBPolyhedron::edge_exists_in_tri(FB_Triangle& tri, int v0, int v1)
{
FB_Edge *edge;
std::vector<FB_Edge*>::iterator dp;

  dp = tri.edge_list.begin();
  while ( dp != tri.edge_list.end() ) {
    edge = *dp;
    if ( ( (edge->v0) == v0 ) && ( (edge->v1) == v1 ) ||
         ( (edge->v0) == v1 ) && ( (edge->v1) == v0 ) )
      return true;
    dp++;  
  }
  
  return false;
}

void FBPolyhedron::add_new_triangle_data()
{
unsigned int i;
FB_Triangle *tri;

  for ( i = original_numtris; i < tris.size(); i++ ) {
  
    tri = tris[i];
//  make the bounding box
    make_tri_boundingbox(tri);
//  make the plane coefficients
    make_tri_plane_coeffs(tri);
  }
  
}

void FBPolyhedron::make_tri_plane_coeffs(FB_Triangle *tri)
{
FB_Coord *mycoord;
double x1, x2, x3, y1, y2, y3, z1, z2, z3, e1x, e1y, e1z, e2x, e2y, e2z;
double a, b, c, d, dtemp;

     mycoord = verts[tri->v0];
     x1 = mycoord->coord[0];
     y1 = mycoord->coord[1];
     z1 = mycoord->coord[2];
     mycoord = verts[tri->v1];
     x2 = mycoord->coord[0];
     y2 = mycoord->coord[1];
     z2 = mycoord->coord[2];
     mycoord = verts[tri->v2];
     x3 = mycoord->coord[0];
     y3 = mycoord->coord[1];
     z3 = mycoord->coord[2];
     e1x = x1 - x2; e1y = y1 - y2; e1z = z1 - z2;
     e2x = x3 - x2; e2y = y3 - y2; e2z = z3 - z2;
     a = e1z*e2y - e2z*e1y;
     b = e1x*e2z - e2x*e1z;
     c = e1y*e2x - e2y*e1x;
     dtemp = sqrt(a*a + b*b + c*c);
     if ( dtemp > EPSILON2 ) {
       a /= dtemp;
       b /= dtemp;
       c /= dtemp;
     } else {
       PRINT_WARNING("small-area triangle\n");
     }
     d = -(a*x1 + b*y1 + c*z1);
     tri->a = a; tri->b = b; tri->c = c; tri->d = d;

}

void FBPolyhedron::make_tri_boundingbox(FB_Triangle *tri)
{
double xmin, ymin, zmin, xmax, ymax, zmax;
int j;
int connections[3];
FB_Coord *mycoord;

     xmin = ymin = zmin = CUBIT_DBL_MAX;
     xmax = ymax = zmax = -xmin;
     connections[0] = tri->v0; connections[1] = tri->v1; connections[2] = tri->v2;
     for ( j = 0; j < 3; j++ ) { // make the bounding box
       mycoord = verts[connections[j]];
       xmin = ( xmin < mycoord->coord[0] ) ? xmin : mycoord->coord[0];
       xmax = ( xmax > mycoord->coord[0] ) ? xmax : mycoord->coord[0];
       ymin = ( ymin < mycoord->coord[1] ) ? ymin : mycoord->coord[1];
       ymax = ( ymax > mycoord->coord[1] ) ? ymax : mycoord->coord[1];
       zmin = ( zmin < mycoord->coord[2] ) ? zmin : mycoord->coord[2];
       zmax = ( zmax > mycoord->coord[2] ) ? zmax : mycoord->coord[2];
     }
     if ( (xmax - xmin) < BOX_CRACK ) {
       xmax += BOX_CRACK; xmin -= BOX_CRACK;
     }
     if ( (ymax - ymin) < BOX_CRACK ) {
       ymax += BOX_CRACK; ymin -= BOX_CRACK;
     }
     if ( (zmax - zmin) < BOX_CRACK ) {
       zmax += BOX_CRACK; zmin -= BOX_CRACK;
     }     
     tri->boundingbox.xmin = xmin; tri->boundingbox.xmax = xmax;
     tri->boundingbox.ymin = ymin; tri->boundingbox.ymax = ymax;
     tri->boundingbox.zmin = zmin; tri->boundingbox.zmax = zmax;

     //  Update the object's bounding box
     polyxmin = ( polyxmin < xmin ) ? polyxmin : xmin;
     polyymin = ( polyymin < ymin ) ? polyymin : ymin;
     polyzmin = ( polyzmin < zmin ) ? polyzmin : zmin;
     polyxmax = ( polyxmax > xmax ) ? polyxmax : xmax;
     polyymax = ( polyymax > ymax ) ? polyymax : ymax;
     polyzmax = ( polyzmax > zmax ) ? polyzmax : zmax;

}

void FBPolyhedron::removeduddedtriangles()  //  Compacts the tris vector
{
unsigned int i;
int j;

  j = 0;
  for ( i = 0; i < tris.size(); i++ ) {
    if ( tris[i]->dudded == true ) continue;
    tris[j] = tris[i]; 
     goodtris.push_back(j);  
    j++;
  }
  tris.resize(j);
}

  //find the largest and smallest angles in this triangle
bool FBPolyhedron::min_max_angles_in_fb_triangle(FB_Triangle *triangle,
                                                double& min_angle,
                                                double& max_angle)
{
  CubitVector vert_0(verts[triangle->v0]->coord[0],
                     verts[triangle->v0]->coord[1],
                     verts[triangle->v0]->coord[2]);
  CubitVector vert_1(verts[triangle->v1]->coord[0],
                     verts[triangle->v1]->coord[1],
                     verts[triangle->v1]->coord[2]);
  CubitVector vert_2(verts[triangle->v2]->coord[0],
                     verts[triangle->v2]->coord[1],
                     verts[triangle->v2]->coord[2]);
  CubitVector sides[3];
  sides[0] = vert_1 - vert_0;
  sides[1] = vert_2 - vert_1;
  sides[2]= vert_0 - vert_2;
  if(sides[0].length_squared() < EPSILON ||
     sides[1].length_squared() < EPSILON ||
     sides[2].length_squared() < EPSILON ){
    min_angle =0.0;
    max_angle =180.0;
    return false;
  }
  double curr_angle;
  min_angle = sides[1].interior_angle(-sides[0]);
  max_angle = min_angle;
  curr_angle = sides[2].interior_angle(-sides[1]);
  if(curr_angle<min_angle){
    min_angle=curr_angle;
  }
  if(curr_angle>max_angle){
    min_angle=curr_angle;
  }
  curr_angle = sides[0].interior_angle(-sides[2]);
  if(curr_angle<min_angle){
    min_angle=curr_angle;
  }
  if(curr_angle>max_angle){
    min_angle=curr_angle;
  }
  return true;
  
}
//draw the "boundary edges" of the polyhedron.
void FBPolyhedron::debug_draw_boundary_edges(int color)
{
  unsigned int i;
  FB_Triangle* triangle=NULL;
  for(i = 0; i<tris.size(); i++){
    if(!tris[i]->dudded){
      triangle = tris[i];
      unsigned int counter = 0;
      CubitVector vert_0(verts[triangle->v0]->coord[0],
                           verts[triangle->v0]->coord[1],
                         verts[triangle->v0]->coord[2]);
      CubitVector vert_1(verts[triangle->v1]->coord[0],
                         verts[triangle->v1]->coord[1],
                         verts[triangle->v1]->coord[2]);
      CubitVector vert_2(verts[triangle->v2]->coord[0],
                         verts[triangle->v2]->coord[1],
                         verts[triangle->v2]->coord[2]);
      if(triangle->cubitedge0index){
        counter++;
        GfxDebug::draw_line(vert_0, vert_1, color);
      }
      if(triangle->cubitedge1index){
        counter++;
        GfxDebug::draw_line(vert_1, vert_2, color);
      }
      if(triangle->cubitedge2index){
        counter++;
        GfxDebug::draw_line(vert_2, vert_0, color);
      }
      if(counter != triangle->edge_list.size()){
        PRINT_WARNING("Possible debug problem.\n");
      }
    }
  }
}

//draw a single triangle
void FBPolyhedron::debug_draw_fb_triangle(FB_Triangle *triangle)
{
  CubitVector vert_0(verts[triangle->v0]->coord[0],
                     verts[triangle->v0]->coord[1],
                     verts[triangle->v0]->coord[2]);
  CubitVector vert_1(verts[triangle->v1]->coord[0],
                     verts[triangle->v1]->coord[1],
                     verts[triangle->v1]->coord[2]);
  CubitVector vert_2(verts[triangle->v2]->coord[0],
                     verts[triangle->v2]->coord[1],
                     verts[triangle->v2]->coord[2]);
  GfxDebug::draw_point(vert_0, CUBIT_RED);
  GfxDebug::draw_point(vert_1, CUBIT_RED);
  GfxDebug::draw_point(vert_2, CUBIT_RED);
  GfxDebug::draw_line(vert_0, vert_1, CUBIT_BLUE);
  GfxDebug::draw_line(vert_1, vert_2, CUBIT_BLUE);
  GfxDebug::draw_line(vert_2, vert_0, CUBIT_BLUE);
  int draw_normal = 1;
  if(draw_normal){
    CubitVector center_pos = (vert_0+vert_1+vert_2)/3.0;
    CubitVector opp_pos = center_pos + (vert_1-vert_0)*(vert_2-vert_0);
    GfxDebug::draw_line(center_pos, opp_pos, CUBIT_WHITE);
  }
}

//get the largest and smallest angles in the polyhedron.  That is,
//find the largest angle in any triangle in the polyhedron and find
//the smallest angle in any triangle in the polyhedron.
bool FBPolyhedron::min_max_angles_in_polyhedron(double& min_angle,
                                                  double& max_angle)
{
  double curr_min;
  double curr_max;
  if(tris.size() < 1){
    min_angle = 0.0;
    max_angle = 180.0;
    return false;
  }
  min_max_angles_in_fb_triangle(tris[0], min_angle, max_angle);
  unsigned int i;
  for(i = 1; i<tris.size(); i++){
    if(!min_max_angles_in_fb_triangle(tris[i], curr_min, curr_max)){
      return false;
    }
    if(curr_min<min_angle){
      min_angle=curr_min;
    }
    if(curr_max>max_angle){
      max_angle=curr_max;
    }
  }
  return true;
}

#ifdef KEEP_BOYD10_KEEP
//function that attempts to remove small angles be flipping edges...
bool FBPolyhedron::remove_small_angles_in_triangle_range(int lower_index,
                                                         int upper_index)
{
    //we currently only handle obtuse, degenerate triangles...
    //see if the triangle falls in this range
  int mydebug = 0;
  const double min_angle_in_degrees = 1;
  const double max_angle_in_degrees = 120;
  int j;
  for(j=lower_index;j<upper_index;j++){
      //this triangle has an extremely small angl
    double min_angle, max_angle;
    if(!min_max_angles_in_fb_triangle(tris[j], min_angle, max_angle))
      return false;
    
    if(min_angle < min_angle_in_degrees && max_angle>max_angle_in_degrees)
    {
        //PRINT_INFO("Found a degenerate, obtuse triangle.\n");
        //two vertices on longest edge of degenerate triangle
      int longest_edge_v1 = 0;
      int longest_edge_v2 = 0;
        //all vertices on triangle
      int v_indices[3];
      v_indices[0]=tris[j]->v0;
      v_indices[1]=tris[j]->v1;
      v_indices[2]=tris[j]->v2;
      double longest_length = -1.0;
      double temp_length = 0.0;
      int k=0;
        //find the longest edge
      for(k=0;k<3;k++){
        CubitVector v1(verts[v_indices[k]]->coord[0],
                       verts[v_indices[k]]->coord[1],
                       verts[v_indices[k]]->coord[2]);
        CubitVector v2(verts[v_indices[(k+1)%3]]->coord[0],
                       verts[v_indices[(k+1)%3]]->coord[1],
                       verts[v_indices[(k+1)%3]]->coord[2]);
        temp_length = (v1-v2).length_squared();
              
        if(temp_length>longest_length){
          longest_length = temp_length;
          longest_edge_v1 = v_indices[k];
          longest_edge_v2 = v_indices[(k+1)%3];
        }
      }
        //look for the other triangle in this range that shares this edge
      int found_it = 0;
      k=0;
      while (!found_it && k < (int)tris.size()){
        if(k!=j && tris[k]->are_fb_coords_in_fb_triangle(
          longest_edge_v1,longest_edge_v2)){
          found_it=1;
        }
        if(!found_it){
          k++;
        }
      }
        //we have found the two triangles adjacent to the longest
        //edge of the degenerate tri.
      if(!found_it){
          //PRINT_WARNING("Didn't find the other triangle\n");
        return false;
      }
//       PRINT_INFO("Found it... j = %i, k = %i\n",j,k);
//       PRINT_INFO("edge lists sizes %i, %i\n",tris[j]->edge_list.size(),tris[k]->edge_list.size());
      double min_angle_before_swap;
      min_max_angles_in_fb_triangle(tris[j], min_angle_before_swap,
                                    max_angle);
      // PRINT_INFO("Inititial tri j angle = %f\n",min_angle);
      min_max_angles_in_fb_triangle(tris[k], min_angle, max_angle);
      if(min_angle<min_angle_before_swap){
        min_angle_before_swap = min_angle;
      }
      if(mydebug){
        PRINT_INFO("Inititial tri k angle = %f\n",min_angle);
        GfxDebug::clear();
        debug_draw_fb_triangle(tris[j]);
        debug_draw_fb_triangle(tris[k]);    
        GfxDebug::mouse_xforms();
      }
        //now do a swap to remove the longest edge...
      int other_vj;
      int other_vk;
      int index_for_j = 0;
      int index_for_k =0;
      int old_value_for_j=0;
      int old_value_for_k=0;
      if(tris[j]->v0 == longest_edge_v1 ||
         tris[j]->v0 == longest_edge_v2 ){
        if(tris[j]->v1 == longest_edge_v1 ||
           tris[j]->v1 == longest_edge_v2){
          other_vj =tris[j]->v2;
          index_for_j=1;
        }
        else{
          other_vj =tris[j]->v1;
          index_for_j=0;
        }
      }
      else{
        other_vj =tris[j]->v0;
        index_for_j=2;
      }
      if(tris[k]->v0 == longest_edge_v1 ||
         tris[k]->v0 == longest_edge_v2 ){
        if(tris[k]->v1 == longest_edge_v1 ||
           tris[k]->v1 == longest_edge_v2){
          other_vk = tris[k]->v2;
          old_value_for_k=tris[k]->v1;
          tris[k]->v1=other_vj;
          index_for_k = 1;
        }
        else{
          other_vk = tris[k]->v1;
          old_value_for_k=tris[k]->v0;
          tris[k]->v0=other_vj;
          index_for_k=0;
        }
      }
      else{
        other_vk = tris[k]->v0;
        old_value_for_k=tris[k]->v2;
        tris[k]->v2=other_vj;
        index_for_k=2;
      }
      if(index_for_j == 0){
        old_value_for_j=tris[j]->v0;
        tris[j]->v0=other_vk;
      }
      else if(index_for_j == 1){
        old_value_for_j=tris[j]->v1;
        tris[j]->v1=other_vk;
      }
      else{
        old_value_for_j=tris[j]->v2;
        tris[j]->v2=other_vk;
      }
      double min_angle_after_swap;
      min_max_angles_in_fb_triangle(tris[j], min_angle_after_swap,
                                    max_angle);
        // PRINT_INFO("Modified tri j angle = %f\n",min_angle);
      min_max_angles_in_fb_triangle(tris[k], min_angle, max_angle);
      if(min_angle<min_angle_after_swap){
        min_angle_after_swap=min_angle;
      }
      if(min_angle_after_swap < min_angle_before_swap){
          // PRINT_WARNING("Swap was unsuccessful in removing bad angle.");
        if(index_for_j == 0){
          tris[j]->v0=old_value_for_j;
        }
        else if(index_for_j == 1){
          
          tris[j]->v1=old_value_for_j;
        }
        else{
          tris[j]->v2=old_value_for_j;
        }
        if(index_for_k == 0){
          tris[k]->v0=old_value_for_k;
        }
        else if(index_for_k == 1){
          
          tris[k]->v1=old_value_for_k;
        }
        else{
          tris[k]->v2=old_value_for_k;
        }
      }
      if(mydebug){
        PRINT_INFO("Modified tri k angle = %f\n",min_angle);
        GfxDebug::clear();
        debug_draw_fb_triangle(tris[j]);
        debug_draw_fb_triangle(tris[k]);    
        GfxDebug::mouse_xforms();
      }
    }
  }
  return true;
  
}
#endif

  
