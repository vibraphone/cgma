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

#include "CubitMessage.hpp"
#include "FBClassify.hpp"
#include "IntegerHash.hpp"
#include "GfxDebug.hpp"
#include <stack>
//make this the same CUBIT_RESABS???
//const double EPSILON_CLASSIFY = 1.e-12;
const double EPSILON_CLASSIFY = 1.e-10;
FBClassify::FBClassify()
{
  number_of_groups = 0;
  polya = polyb = 0;
}

FBClassify::~FBClassify()
{

}

void FBClassify::SetPoly(FBPolyhedron *poly1, FBPolyhedron *poly2)
{
  polya = poly1;
  polyb = poly2;
}

CubitStatus FBClassify::Group(int which)
{
unsigned int i;
int hashvalue, v0, v1, v2, tv0, tv1, tv2;
//const int numhashbins = 101;
int *hasharrayptr, hasharraysize, j, itri;
FBPolyhedron *poly;
const int primes[] = { 307, 601, 1009, 3001, 6007, 10007, 30011, 60013, 100003 };
int numhashbins;

  if ( which == 1 ) poly = polya;
  else poly = polyb;
  
  if ( poly == 0 ) {
    PRINT_ERROR("ERROR:  Polyhedral object undefined in FBClassify::Group\n");
    return CUBIT_FAILURE;
  }

  i = poly->verts.size();
  if ( i < 3000 ) numhashbins = primes[0];
  else if ( i < 6000 ) numhashbins = primes[1];
  else if ( i < 10000 ) numhashbins = primes[2];
  else if ( i < 30000 ) numhashbins = primes[3];
  else if ( i < 60000 ) numhashbins = primes[4];
  else if ( i < 100000 ) numhashbins = primes[5];
  else if ( i < 300000 ) numhashbins = primes[6];
//  else if ( i < 600000 ) numhashbins = primes[7];
  else numhashbins = primes[7]; 

  IntegerHash *hash = new IntegerHash(numhashbins,50);
  e0 = new int[poly->tris.size()];
  e1 = new int[poly->tris.size()];
  e2 = new int[poly->tris.size()];

  for ( i = 0; i < poly->tris.size(); i++ ) {
    e0[i] = e1[i] = e2[i]= NO_EDGE_NBR; // initialize 
    group.push_back(UNSET);
    v0 = poly->tris[i]->v0;
    v1 = poly->tris[i]->v1;
    v2 = poly->tris[i]->v2;
    //  Put the triangle sequence, i, in the hash list for each of the 3 edges.
    hash->addtoHashList((v0+v1)%numhashbins,i);
    hash->addtoHashList((v1+v2)%numhashbins,i);
    hash->addtoHashList((v2+v0)%numhashbins,i);
  }
/*  int jmin, jmax, jave;
  jmin = 100000000;
  jmax = -100000000;
  jave = 0;
  for ( i = 0; i < poly->tris.size(); i++ ) {
    hasharrayptr = hash->getHashBin(i,&hasharraysize);
    if ( hasharraysize < jmin ) jmin = hasharraysize;
    if ( hasharraysize > jmax ) jmax = hasharraysize;
    jave += hasharraysize;
    
  } 
  jave /= poly->tris.size();
  printf("jmin jmax jave %d %d %d\n",jmin, jmax, jave);
  */
  
  for ( i = 0; i < poly->tris.size(); i++ ) {
    v0 = poly->tris[i]->v0;
    v1 = poly->tris[i]->v1;
    v2 = poly->tris[i]->v2;
    hashvalue = (v0+v1)%numhashbins;  // get the hash value for edge 0
    hasharrayptr = hash->getHashBin(hashvalue,&hasharraysize);
    for ( j = 0; j < hasharraysize; j++ ) {  
      //  Go through the list and find the other triangle, itri, for each edge.
      //  Then assign its sequence to the edge array.
      itri = hasharrayptr[j];
      if ( (unsigned int)itri == i ) continue;
      tv0 = poly->tris[itri]->v0;
      tv1 = poly->tris[itri]->v1;
      tv2 = poly->tris[itri]->v2;
      if ( ((v1 == tv0) && (v0 == tv1)) || ((v0 == tv0) && (v1 == tv1)) ) {
        e0[i] = itri;
      } else if ( ((v1 == tv1) && (v0 == tv2)) || ((v0 == tv1) && (v1 == tv2)) ) {
        e0[i] = itri;
      } else if ( ((v1 == tv2) && (v0 == tv0)) || ((v0 == tv2) && (v1 == tv0)) ) {
        e0[i] = itri;
      } 
    } 
    hashvalue = (v1+v2)%numhashbins;
    hasharrayptr = hash->getHashBin(hashvalue,&hasharraysize);
    for ( j = 0; j < hasharraysize; j++ ) {
      itri = hasharrayptr[j];
      if ( (unsigned int)itri == i ) continue;
      tv0 = poly->tris[itri]->v0;
      tv1 = poly->tris[itri]->v1;
      tv2 = poly->tris[itri]->v2;
      if ( ((v1 == tv0) && (v2 == tv1)) || ((v2 == tv0) && (v1 == tv1)) ) {
        e1[i] = itri;
      } else if ( ((v1 == tv1) && (v2 == tv2)) || ((v2 == tv1) && (v1 == tv2)) ) {
        e1[i] = itri;
      } else if ( ((v1 == tv2) && (v2 == tv0)) || ((v2 == tv2) && (v1 == tv0)) ) {
        e1[i] = itri;
      } 
    } 
    hashvalue = (v2+v0)%numhashbins;
    hasharrayptr = hash->getHashBin(hashvalue,&hasharraysize);
    for ( j = 0; j < hasharraysize; j++ ) {
      itri = hasharrayptr[j];
      if ( (unsigned int)itri == i ) continue;
      tv0 = poly->tris[itri]->v0;
      tv1 = poly->tris[itri]->v1;
      tv2 = poly->tris[itri]->v2;
      if ( ((v0 == tv0) && (v2 == tv1)) || ((v2 == tv0) && (v0 == tv1)) ) {
        e2[i] = itri;
      } else if ( ((v0 == tv1) && (v2 == tv2)) || ((v2 == tv1) && (v0 == tv2)) ) {
        e2[i] = itri;
      } else if ( ((v0 == tv2) && (v2 == tv0)) || ((v2 == tv2) && (v0 == tv0)) ) {
        e2[i] = itri;
      }
    }       
  }
  //  Now we have to remove the other-side triangles where there was an intersection
  //  edge. 
  for ( i = 0; i < poly->intersection_edges.size(); i++ ) {
    v0 = poly->intersection_edges[i]->v0;
    v1 = poly->intersection_edges[i]->v1;
    hashvalue = (v0+v1)%numhashbins;
    hasharrayptr = hash->getHashBin(hashvalue,&hasharraysize);
    for ( j = 0; j < hasharraysize; j++ ) {
      itri = hasharrayptr[j];
      tv0 = poly->tris[itri]->v0;
      tv1 = poly->tris[itri]->v1;
      tv2 = poly->tris[itri]->v2;
      if ( ((v0 == tv0) && (v1 == tv1)) || ((v0 == tv1) && (v1 == tv0)) ) {
        e0[itri] = NO_EDGE_NBR;
      };
      if ( ((v0 == tv0) && (v1 == tv2)) || ((v0 == tv2) && (v1 == tv0)) ) {
        e2[itri] = NO_EDGE_NBR;
      };
      if ( ((v0 == tv1) && (v1 == tv2)) || ((v0 == tv2) && (v1 == tv1)) ) {
        e1[itri] = NO_EDGE_NBR;
      }
    }    
  }
    
  //  Group the triangles that are neighbors.
  for ( i = 0; i < poly->tris.size(); i++ ) {
    if ( group[i] == UNSET ) {
      fill_group(i,number_of_groups++);
    }
  }

  delete[] e0; delete[] e1; delete[] e2; 
  delete hash;
    
  return CUBIT_SUCCESS;
}

void FBClassify::fill_group(int itri, int ngroup)
{
std::stack<int> vstack;
int ktri;

  group[itri] = ngroup;
  
  if ( (e0[itri] != NO_EDGE_NBR) && (group[e0[itri]] == UNSET) ) 
    vstack.push(e0[itri]); 
  if ( (e1[itri] != NO_EDGE_NBR) && (group[e1[itri]] == UNSET) ) 
    vstack.push(e1[itri]); 
  if ( (e2[itri] != NO_EDGE_NBR) && (group[e2[itri]] == UNSET) ) 
    vstack.push(e2[itri]); 
  while (vstack.size() > 0 ) {
    ktri = vstack.top();
    vstack.pop();
    group[ktri] = ngroup;
    if ( (e0[ktri] != NO_EDGE_NBR) && (group[e0[ktri]] == UNSET) ) 
      vstack.push(e0[ktri]); 
    if ( (e1[ktri] != NO_EDGE_NBR) && (group[e1[ktri]] == UNSET) ) 
      vstack.push(e1[ktri]); 
    if ( (e2[ktri] != NO_EDGE_NBR) && (group[e2[ktri]] == UNSET) ) 
      vstack.push(e2[ktri]); 
  }
}

CubitStatus FBClassify::CharacterizeGroups(int which, bool other_is_planar)
{
  int numberdone;
  unsigned int i;
  bool ifoundit;
  int itri = -1, type;
  FBPolyhedron *poly;

  if ( which == 1 ) poly = polya;
  else poly = polyb;
  
  if ( poly == 0 ) {
    PRINT_ERROR("ERROR:  Polyhedral object undefined in FBClassify::Group\n");
    return CUBIT_FAILURE;
  }
  
  numberdone = 0;
  i = 0; 
  while ( numberdone < number_of_groups ) {
    ifoundit = false;
    while ( i < poly->tris.size() ) {
      if ( group[i] == numberdone ) {
        ifoundit = true;
        itri = (int)i;
        break;      
      }
      i++;
    }
    if ( ifoundit == true ) {
      if ( other_is_planar == false ) 
        type = classify(itri, which);
      else 
        type = classify_against_plane(itri, which);      
      group_characterization.push_back(type);
      numberdone++;
    } else {
      PRINT_ERROR("Error in FBClassify::CharacterizeGroups\n");
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

int FBClassify::classify_against_plane(int itri, int which)
{
FBPolyhedron *polyref, *polyobj; 
int type, v0, v1, v2;
double xbary, ybary, zbary, a, b, c;
;

  type = FB_ORIENTATION_UNDEFINED;
  if ( which == 1 ) {
    polyobj = polyb;
    polyref = polya;
  } else if ( which == 2 ) {
    polyobj = polya;
    polyref = polya;
  } else {
    PRINT_ERROR("ERROR in FBClassify::classify\n");
    return type;
  }
  v0 = polyref->tris[itri]->v0;
  v1 = polyref->tris[itri]->v1;
  v2 = polyref->tris[itri]->v2;

  xbary = ( polyref->verts[v0]->coord[0] + 
            polyref->verts[v1]->coord[0] + 
            polyref->verts[v2]->coord[0] )/3.;
  ybary = ( polyref->verts[v0]->coord[1] + 
            polyref->verts[v1]->coord[1] + 
            polyref->verts[v2]->coord[1] )/3.;
  zbary = ( polyref->verts[v0]->coord[2] + 
            polyref->verts[v1]->coord[2] + 
            polyref->verts[v2]->coord[2] )/3.;
  a = polyref->tris[itri]->a;
  b = polyref->tris[itri]->b;
  c = polyref->tris[itri]->c;

  //  Figure out which side of the plane we are on.  Since all
  //  of the plane's triangles have the same plane equation
  //  coefficients, might as well use the first one.
  
double obj_tri_a, obj_tri_b, obj_tri_c, obj_tri_d, dotprod, disttoplane;

  obj_tri_a = polyobj->tris[0]->a;
  obj_tri_b = polyobj->tris[0]->b;
  obj_tri_c = polyobj->tris[0]->c;
  obj_tri_d = polyobj->tris[0]->d;
      
  disttoplane = obj_tri_a*xbary + obj_tri_b*ybary + obj_tri_c*zbary + obj_tri_d;    
  if ( disttoplane > EPSILON ) return FB_ORIENTATION_OUTSIDE;
  else if ( disttoplane < -EPSILON ) return FB_ORIENTATION_INSIDE; 
  
  dotprod = obj_tri_a*a + obj_tri_b*b + obj_tri_c*c;
  if ( dotprod > 0. ) return FB_ORIENTATION_SAME;
  else return FB_ORIENTATION_OPPOSITE;
  
  return type;
}

//returns an orientation for the triangle relative to the other body.
//This triangle can be inside or outside the other body.  Or, it can
//be opposite or same.  Opposite means this triangle is very close to
// a triangle in the other body and it has an opposite pointing normal.
// Same means it is very close to a trianglein the other body and it
// their normals point in the same direction.
int FBClassify::classify(int itri, int which)
{
    //  "which" is the object, either 1 or 2, that the triangle itri
    // belongs to.
    //  Th e object to test for inside or outside is the other object.
  FBPolyhedron *polyref, *polyobj; 
    //  polyref = itri's polyhedron; polyobj = the other object
  int type, v0, v1, v2;
  double xbary, ybary, zbary, a, b, c;

  type = FB_ORIENTATION_UNDEFINED;
  if ( which == 1 ) {
    polyobj = polyb;
    polyref = polya;
  } else if ( which == 2 ) {
    polyobj = polya;
    polyref = polyb;
  } else {
    PRINT_ERROR("ERROR in FBClassify::classify\n");
    return type;
  }
  int mydebug = 0;
  if(mydebug)
    polyref->debug_draw_fb_triangle(polyref->tris[itri]);

  v0 = polyref->tris[itri]->v0;
  v1 = polyref->tris[itri]->v1;
  v2 = polyref->tris[itri]->v2;

  xbary = ( polyref->verts[v0]->coord[0] + 
            polyref->verts[v1]->coord[0] + 
            polyref->verts[v2]->coord[0] )/3.;
  ybary = ( polyref->verts[v0]->coord[1] + 
            polyref->verts[v1]->coord[1] + 
            polyref->verts[v2]->coord[1] )/3.;
  zbary = ( polyref->verts[v0]->coord[2] + 
            polyref->verts[v1]->coord[2] + 
            polyref->verts[v2]->coord[2] )/3.;
  a = polyref->tris[itri]->a;
  b = polyref->tris[itri]->b;
  c = polyref->tris[itri]->c;

  unsigned int i, num_perturb;
  double obj_tri_a, obj_tri_b, obj_tri_c, obj_tri_d, dotprod;
  double distance_to_other_sqr, closest_distance_to_plane, t;
  double closest_distance_to_other_sqr;
  double xint, yint, zint, distance_to_plane, closest_dotproduct;
  bool perturb, done, foundone;
  double other_xbar, other_ybar, other_zbar;
  int other_tri_0, other_tri_1, other_tri_2;

  perturb = false;
  num_perturb = 0;
  done = false;
  
  while ( (done == false) && (num_perturb < 20) ) {
    closest_dotproduct = -CUBIT_DBL_MAX + 1.;
    closest_distance_to_plane = CUBIT_DBL_MAX;
    closest_distance_to_other_sqr = CUBIT_DBL_MAX;
    foundone = false;
    for ( i = 0; i < polyobj->tris.size(); i++ ) {
      obj_tri_a = polyobj->tris[i]->a;
      obj_tri_b = polyobj->tris[i]->b;
      obj_tri_c = polyobj->tris[i]->c;
      obj_tri_d = polyobj->tris[i]->d;
  
      dotprod = obj_tri_a*a + obj_tri_b*b + obj_tri_c*c;
        //calculate the distance to the other triangles plane
      distance_to_plane = (obj_tri_a*xbary + obj_tri_b*ybary +
                           obj_tri_c*zbary + obj_tri_d);
       
      
      if ( fabs(dotprod) < EPSILON_CLASSIFY ) {
          //  Is the point in the plane?
        if ( fabs(distance_to_plane) < EPSILON_CLASSIFY ) {
            //  Perturb the ray and recast.
          perturb = true;
          num_perturb += 1;
          break;
        }
        continue;
      }
      
      t =-(distance_to_plane)/dotprod;
      if ( t < -EPSILON_CLASSIFY ) continue;
      xint = xbary + a*t;
      yint = ybary + b*t;
      zint = zbary + c*t;

        //  Check whether the intersection point lies in or on
        //  the object triangle's
        //  bounding box.
      if ( (polyobj->tris[i]->boundingbox.xmin - EPSILON > xint) || 
           (polyobj->tris[i]->boundingbox.xmax + EPSILON < xint) ||
           (polyobj->tris[i]->boundingbox.ymin - EPSILON > yint) || 
           (polyobj->tris[i]->boundingbox.ymax + EPSILON < yint) ||
           (polyobj->tris[i]->boundingbox.zmin - EPSILON > zint) || 
           (polyobj->tris[i]->boundingbox.zmax + EPSILON < zint) ) 
        continue;
      
        //  Is the point (xint, yint, zint) inside or on the triangle?
        //  Get a principal projection to make this a 2D problem.
      double xp1, yp1, xp2, yp2, xp3, yp3, ptx, pty;
      int retval;

      if ( (fabs(obj_tri_b) >= fabs(obj_tri_a)) && 
           (fabs(obj_tri_b) >= fabs(obj_tri_c)) ) {
        xp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[0];
        yp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[2];
        xp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[0];
        yp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[2];
        xp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[0];
        yp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[2];
        ptx = xint;
        pty = zint;        
      } else if ( fabs(obj_tri_a) >= fabs(obj_tri_c) ) {
        xp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[1];
        yp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[2];
        xp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[1];
        yp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[2];
        xp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[1];
        yp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[2];
        ptx = yint;
        pty = zint;   
      } else {
        xp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[0];
        yp1 = polyobj->verts[polyobj->tris[i]->v0]->coord[1];
        xp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[0];
        yp2 = polyobj->verts[polyobj->tris[i]->v1]->coord[1];
        xp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[0];
        yp3 = polyobj->verts[polyobj->tris[i]->v2]->coord[1];
        ptx = xint;
        pty = yint;  
      }
      retval = pt_in_tri_2d(ptx,pty,xp1,yp1,xp2,yp2,xp3,yp3);
      if ( (retval == FB_ORIENTATION_INSIDE) ||
           (retval == FB_ORIENTATION_ON) ) {
          //calculate the distance to the other triangle's centroid
        other_tri_0 = polyobj->tris[i]->v0;
        other_tri_1 = polyobj->tris[i]->v1;
        other_tri_2 = polyobj->tris[i]->v2;
        other_xbar = ( polyobj->verts[other_tri_0]->coord[0] + 
                       polyobj->verts[other_tri_1]->coord[0] + 
                       polyobj->verts[other_tri_2]->coord[0] )/3.;
        other_ybar = ( polyobj->verts[other_tri_0]->coord[1] + 
                       polyobj->verts[other_tri_1]->coord[1] + 
                       polyobj->verts[other_tri_2]->coord[1] )/3.;
        other_zbar = ( polyobj->verts[other_tri_0]->coord[2] + 
                       polyobj->verts[other_tri_1]->coord[2] + 
                       polyobj->verts[other_tri_2]->coord[2] )/3.;
        
          //calculate the distance (squared) to the other triangle's centroid
        distance_to_other_sqr = ( (xbary-other_xbar)*(xbary-other_xbar) +
                                  (ybary-other_ybar)*(ybary-other_ybar) +
                                  (zbary-other_zbar)*(zbary-other_zbar) );
          //if this is the closest other triangle so far...
        if(closest_distance_to_other_sqr > distance_to_other_sqr){
            //then we found one, and update the closest distance, dot prod,
            // and distance to other plane.
          foundone = true;
          closest_distance_to_other_sqr = distance_to_other_sqr;
          closest_dotproduct = dotprod;
          if(mydebug){
            polyobj->debug_draw_fb_triangle(polyobj->tris[i]);
            GfxDebug::mouse_xforms();
          }
          closest_distance_to_plane = distance_to_plane;
            //  This is the closest triangle.
          if ( fabs(closest_distance_to_plane) < EPSILON_CLASSIFY ) 
            break;   
        }
      }       
    }
    if ( perturb == false ) done = true;
    else {
        //  perturb the ray and try again.
      perturb_the_ray(xbary, ybary, zbary);      
      perturb = false;
    }
  }
    //if we are very close to the plane of the closest other triangle
    // then we are either going to classify as opposite or same depending
    // on the relationship of the normals
  if ( (fabs(closest_distance_to_plane) < EPSILON_CLASSIFY ) ){
    if ( closest_dotproduct > 0.0 )
      type = FB_ORIENTATION_SAME;
    else
      type = FB_ORIENTATION_OPPOSITE;
  }
    //otherwise we are going to classify as inside or outside.  If we
    // didn't find any triangles that projected into our triangle,
    // we must be outside.  Otherwise we compare the normals to determine
    // whether we are inside or outside.
  else {
    if  ( foundone == false )
      type = FB_ORIENTATION_OUTSIDE;
    else if ( closest_dotproduct > 0.0 )
      type = FB_ORIENTATION_INSIDE;
    else
      type = FB_ORIENTATION_OUTSIDE;  
  }
  if(mydebug)
    GfxDebug::display_all();
  return type;
}

void FBClassify::perturb_the_ray(double &xbary, double &ybary, double &zbary)
{
  xbary += 1.e-4*(double(rand())/(RAND_MAX+1.0)-0.5);
  ybary += 1.e-4*(double(rand())/(RAND_MAX+1.0)-0.5);
  zbary += 1.e-4*(double(rand())/(RAND_MAX+1.0)-0.5);
}

int FBClassify::pt_in_tri_2d(double xpt, double ypt,
                              double x0, double y0,
		              double x1, double y1,
		              double x2, double y2)
{
//  From Schneider & Eberly, "Geometric Tools for COmputer Graphics",
//  Chap. 13.3.1.  If triangle is needle-thin, CUBIT_FAILURE might be
//  returned, in wich case is_point_in is undefined.

double c0, c1, c2;
double e0x, e1x, e2x, e0y, e1y, e2y;
double n0x, n1x, n2x, n0y, n1y, n2y;
double denom0, denom1, denom2;
int result;

  e0x = x1 - x0; e0y = y1 - y0;  
  e1x = x2 - x1; e1y = y2 - y1;  
  e2x = x0 - x2; e2y = y0 - y2;  
  n0x = e0y; n0y = -e0x;
  n1x = e1y; n1y = -e1x;
  n2x = e2y; n2y = -e2x;
  denom0 = n1x*e0x + n1y*e0y;
  if ( fabs(denom0) < EPSILON_CLASSIFY ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return FB_ORIENTATION_UNDEFINED;
  }
  denom1 = n2x*e1x + n2y*e1y;
  if ( fabs(denom1) < EPSILON_CLASSIFY ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return FB_ORIENTATION_UNDEFINED;
  }
  denom2 = n0x*e2x + n0y*e2y;
  if ( fabs(denom2) < EPSILON_CLASSIFY ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return FB_ORIENTATION_UNDEFINED;
  }
  
  c0 = -( n1x*(xpt-x1) + n1y*(ypt-y1) )/denom0;
  c1 = -( n2x*(xpt-x2) + n2y*(ypt-y2) )/denom1;
  c2 = -( n0x*(xpt-x0) + n0y*(ypt-y0) )/denom2;

  if ( (c0 > 0.0) && (c1 > 0.0) && (c2 > 0.0) ) result = FB_ORIENTATION_INSIDE;
  else if ( (c0 < 0.0) || (c1 < 0.0) || (c2 < 0.0) ) result = FB_ORIENTATION_OUTSIDE;
  else result = FB_ORIENTATION_ON;

  return result;

}

void FBClassify::get_group(std::vector<int> **this_group,
                           std::vector<int> **this_group_characterization)
{
  *this_group = &group;
  *this_group_characterization = &group_characterization;
}
