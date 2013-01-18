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

#include <math.h>
#include <algorithm>
#include <vector>
#include "CubitMessage.hpp"
#include "FBTiler.hpp"

FBTiler::FBTiler(std::vector<FB_Coord *>& my_verts, int pd, int sd, int sequence,
                 double a, double b, double c,
                 std::vector<int> *tri_list)
{

  verts = my_verts;
  p_dir = pd;
  s_dir = sd;
  xnorm = a;
  ynorm = b;
  znorm = c;
  parent = sequence;
  my_tri_list = tri_list;
}

FBTiler::~FBTiler()
{

}

CubitStatus FBTiler::Tile_Poly(std::vector<int> *coordlist)
{
int v0, v1, v2, v3;
std::vector<int>::iterator it, itmin, itmax;
double min_p, min_s, max_p, max_s, test_p, test_s;
CubitStatus status;

  status = CUBIT_SUCCESS;
  
  it = coordlist->begin();
  
  if ( coordlist->size() == 3 ) {
    v1 = *it++;  v2 = *it++; v3 = *it++;  
    add_triangle(v1, v2, v3);
    return status;
  }

  //  Get location in list of min and max coord in p_dir.
  min_p = min_s = CUBIT_DBL_MAX;
  max_p = max_s = -CUBIT_DBL_MAX;
  while ( it != coordlist->end() ) {
    v0 = *it;
    test_p = verts[v0]->coord[p_dir];
    test_s = verts[v0]->coord[s_dir];
//    if ( test_p < min_p ) {
    if ( (test_p - min_p) <= -EPSILON ) {
      itmin = it;
      min_p = test_p;
      if ( test_s < min_s ) 
        min_s = test_s;
    } 
//    if ( test_p > max_p ) {
    if ( (test_p - max_p) >= EPSILON ) {
      itmax = it;
      max_p = test_p;
      if ( test_s > max_s )
        max_s = test_s;      
    }
    if ( fabs(test_p - min_p) < EPSILON ) {
      if ( test_s < min_s ) {
        itmin = it;
        min_s = test_s;
      }
    } 
    if ( fabs(test_p - max_p) < EPSILON ) {
      if ( test_s > max_s ) {
        itmax = it;
        max_s = test_s;
      }    
    }

    it++;
  }  

  //  Now we've got to start at the min coord and make a single sorted list  
  //  of the coords, from min to max, noting whether the coord is on the
  //  LEFT or RIGHT chain, or -- for the top and bottom coords -- on BOTH chains.

  sortedchainlist.clear();  
  
  std::vector<int>::iterator itleft, itright, itlistbegin, itlistend;
  it = itleft = itright = itmin;  
  itlistbegin = coordlist->begin();
  itlistend = coordlist->end();
  
  decrement_list_ptr(itleft,itlistbegin,itlistend);
  increment_list_ptr(itright,itlistbegin,itlistend); 
  
  v0 = *itmin;
  int v0left, v0right;
  double leftpdircoord, rightpdircoord;
  FBTilerChainVert *cvert;
  cvert = new FBTilerChainVert(v0,BOTH);
  sortedchainlist.push_back(cvert);
    v0left = *itleft; v0right = *itright;  
    leftpdircoord = verts[v0left]->coord[p_dir];
    rightpdircoord = verts[v0right]->coord[p_dir];

  while(1) {  
    while ( (leftpdircoord - rightpdircoord) > -EPSILON2 ) {
      if ( itright == itmax ) break;
      v0 = v0right;
      cvert = new FBTilerChainVert(v0,LEFT);
      sortedchainlist.push_back(cvert);
      increment_list_ptr(itright,itlistbegin,itlistend); 
      v0right = *itright;
      rightpdircoord = verts[v0right]->coord[p_dir];

    }
    while ( (rightpdircoord - leftpdircoord) > -EPSILON2 ) {
      if ( itleft == itmax ) break;
      v0 = v0left;
      cvert = new FBTilerChainVert(v0,RIGHT);
      sortedchainlist.push_back(cvert);
      decrement_list_ptr(itleft,itlistbegin,itlistend);
      v0left = *itleft;
      leftpdircoord = verts[v0left]->coord[p_dir];

    }    
    if ( (itright == itmax) || (itleft == itmax) ) break;
  }

//  In case we have some horizontal edges (wrt p_dir) left on the list,
//  the next two while statements will get them.
  while ( itright != itmax ) {
      v0 = v0right;
      cvert = new FBTilerChainVert(v0,LEFT);
      sortedchainlist.push_back(cvert);
      increment_list_ptr(itright,itlistbegin,itlistend); 
      v0right = *itright;
  } 
  while ( itleft != itmax ) {
       v0 = v0left;
      cvert = new FBTilerChainVert(v0,RIGHT);
      sortedchainlist.push_back(cvert);
      decrement_list_ptr(itleft,itlistbegin,itlistend);
      v0left = *itleft;
  } 

//  Put the top coord on the list.

  v0 = *itmax;
  cvert = new FBTilerChainVert(v0,BOTH);
  sortedchainlist.push_back(cvert);
  status = retriangulate(coordlist);
std::vector<FBTilerChainVert*>::iterator itv, itvlast;  

  itvlast = sortedchainlist.end();
  itv = sortedchainlist.begin();
  while ( itv != itvlast ) {
    delete *itv;   
    itv++;
  }
  sortedchainlist.clear();
  return status;
}

int FBTiler::add_triangle(int v1, int v2, int v3)
{
int status;

  status = 1;

  //  See if the winding is OK by comparing the normal to the parent's normal.
  double xn, yn, zn;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double ux, uy, uz, vx, vy, vz, product;

  x1 = verts[v1]->coord[0]; y1 = verts[v1]->coord[1]; z1 = verts[v1]->coord[2];
  x2 = verts[v2]->coord[0]; y2 = verts[v2]->coord[1]; z2 = verts[v2]->coord[2];
  x3 = verts[v3]->coord[0]; y3 = verts[v3]->coord[1]; z3 = verts[v3]->coord[2];
  ux = x2 - x1; uy = y2 - y1; uz = z2 - z1;
  vx = x3 - x1; vy = y3 - y1; vz = z3 - z1;
  xn = uy*vz - uz*vy; yn = uz*vx - ux*vz; zn = ux*vy - uy*vx;
  product = xnorm*xn + ynorm*yn + znorm*zn;
  if ( product < 0.0 ) {  //  Must reverse the connections.
    int itemp;
    itemp = v1; v1 = v2; v2 = itemp;
  }
  my_tri_list->push_back(v1);
  my_tri_list->push_back(v2);
  my_tri_list->push_back(v3);
 
  return status;

}

void FBTiler::dud_from_coord_list(int val, std::vector<int> *ivec)
{
std::vector<int>::iterator itv;
bool ifoundit;
//  Also decrements ivec->size() if val was found.

  ifoundit = false;
  itv = ivec->begin();
  do {
    if ( ifoundit == false ) {
      if ( *itv == val ) {
        ifoundit = true;
        break;  
      }
    }
    itv++;
  } while ( itv != ivec->end() );
  if ( ifoundit == true ) {
    ivec->erase(itv);
  } 

}

bool FBTiler::reflex_angle(int v0, int v1, int v2, std::vector<int> *coordlist)
{
double v0x, v0y, v1x, v1y, v2x, v2y, xbary, ybary;

  v0x = verts[v0]->coord[s_dir]; 
  v0y = verts[v0]->coord[p_dir]; 
  v1x = verts[v1]->coord[s_dir]; 
  v1y = verts[v1]->coord[p_dir]; 
  v2x = verts[v2]->coord[s_dir]; 
  v2y = verts[v2]->coord[p_dir]; 
//  Are the three points colinear?

  if ( fabs((v1x-v0x)*(v2y-v0y) - (v2x-v0x)*(v1y-v0y)) < EPSILON2 ) 
    return true;
        
//  Get coordinates of barycenter of the triangle formed by 
//  v0, v1, and v2

  xbary = (v0x + v1x + v2x)/3.;
  ybary = (v0y + v1y + v2y)/3.;
  
//  Do a point-in-polygon test on this point and the polygon formed
//  by the chain.  If point is in polygon, angle is not reflex; 
//  otherwise, angle is reflex.  Point in polygon algorithm from
//  Schneider and Eberly, Geometric Tools for COmputer Graphics, 
//  Sec. 13.3.

bool inside = false;
std::vector<int>::iterator itv;
int u0, u1;

  itv = coordlist->begin();
  while ( itv != coordlist->end() ) {
    u0 = *itv++;
    if ( itv == coordlist->end() )
      u1 = *(coordlist->begin());
    else 
      u1 = *itv;
    if ( ybary < verts[u1]->coord[p_dir] ) {
    //  u1 is above the hroizontal ray in s_dir from the barycenter.
      if ( verts[u0]->coord[p_dir] <= ybary ) {
      //  u0 is on or below the ray.
        if ( ( (ybary - verts[u0]->coord[p_dir])*
	       (verts[u1]->coord[s_dir] -verts[ u0]->coord[s_dir]) ) >
	     ( (xbary - verts[u0]->coord[s_dir])*
	       (verts[u1]->coord[p_dir] - verts[u0]->coord[p_dir]) ) )
	  inside = !inside; 
	}    
    } else if ( ybary < verts[u0]->coord[p_dir] ) {
    //  U1 is on or below the ray; u0 is above the ray.
        if ( ( (ybary - verts[u0]->coord[p_dir])*
	       (verts[u1]->coord[s_dir] - verts[u0]->coord[s_dir]) ) <
	     ( (xbary - verts[u0]->coord[s_dir])*
	       (verts[u1]->coord[p_dir] - verts[u0]->coord[p_dir]) ) )
	  inside = !inside; 
    }   
  }

  return !inside;  
  
}

CubitStatus FBTiler::retriangulate(std::vector<int> *coordlist)
{
int stacktopadjacency,  this_adjacency, adjacency_case;
std::vector<FBTilerChainVert*> vstack;
std::vector<FBTilerChainVert*>::iterator itv, itvlast;
int v0, v1, v2;
unsigned int i, numtris, totaltris;

//  Put the first two verts on vstack.
  itv = sortedchainlist.begin();
  vstack.push_back(*itv++);  
  vstack.push_back(*itv);
  stacktopadjacency = (*itv)->whichchain;
  itvlast = itv;
  numtris = 0;
  totaltris = sortedchainlist.size() - 2;
  
  while(1) {
    
    if ( numtris == totaltris ) break;
    *itv++;
    if ( itv == sortedchainlist.end() ) {
      PRINT_ERROR("Tiler error in FacetBoolean:  premature end of list.\n");
      return CUBIT_FAILURE;
    }
 
    this_adjacency = (*itv)->whichchain;
    v0 = (*itv)->v0;

    adjacency_case = get_adjacency(stacktopadjacency,this_adjacency);

    switch ( adjacency_case ) {
  
      case 1:
      {
        v1 = vstack[0]->v0;
        for ( i = 1; i < vstack.size(); i++ ) {
          v2 = vstack[i]->v0;
          add_triangle(v0,v1,v2);
          dud_from_coord_list(v1,coordlist);
          numtris += 1;
          v1 = v2;
        }   
        vstack.clear();
        vstack.push_back(*itvlast);
        vstack.push_back(*itv);
        stacktopadjacency = (*itv)->whichchain;
        itvlast = itv;

        break;
      }

      case 2:
      {
        while ( vstack.size() > 1 ) {
      int size = vstack.size() - 1;
      v1 = vstack[size]->v0;
//      int v1chain = vstack[size]->whichchain;
      v2 = vstack[size-1]->v0;
//      bool isreflex = reflex_angle(v0,v1,v2,v1chain);
      bool isreflex = reflex_angle(v0,v1,v2,coordlist);
      if ( isreflex == true ) break;
      add_triangle(v0,v1,v2);
      dud_from_coord_list(v1,coordlist);
      numtris += 1;
      vstack.pop_back();

   }
    
   vstack.push_back(*itv);
   stacktopadjacency = (*itv)->whichchain;
        itvlast = itv;
    
    
        break;
      }  
  
    }
 
  } // endwhile

  return CUBIT_SUCCESS;
  
}


