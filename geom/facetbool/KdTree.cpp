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

#include "KdTree.hpp"
#include <stack>
#include <vector>
#include <math.h>

const int XCUT = 0;
const int YCUT = 1;
const int ZCUT = 2;

KDTree::KDTree()
{
  epsilonkd = 1.e-6;
}

KDTree::~KDTree()
{
int i;
FSBoundingBox* fsbb;
 for ( i = 0; i < 2*numtris-1; i++ ) {
    fsbb = treebox[i];
    delete fsbb;
  }
  delete [] treebox;
  delete [] nextbranch; 
}

int KDTree::makeKDTree(int npoly, const FSBOXVECTOR& boxlist)
{
double *xcenterdata, *ycenterdata, *zcenterdata, *whichcut[3], *cutarray; 
int i, cuttingdir, min_tri_sequence_value, max_tri_sequence_value, 
    median_value;
double xcen, ycen, zcen;
std::vector<TreeStack* > mytreestack;
int icnt;
TreeStack *thistreestack;
FSBoundingBox *thisbox;

  numtris = npoly;
  boxlistptr = boxlist;
  isequence = new int[npoly];
  xcenterdata = new double[npoly];
  ycenterdata = new double[npoly];
  zcenterdata = new double[npoly];
  nextbranch = new int[2*npoly];
  treebox = new FSBoundingBox*[2*npoly];
  
//  Get the arrays of bounding box center points for x, y, and z.  We will use 
//  these arrays to find the median values along each direction.
  for ( i = 0; i < npoly; i++ ) {
    isequence[i] = i;
    xcen = 0.5*(boxlist[i]->xmin + boxlist[i]->xmax);
    ycen = 0.5*(boxlist[i]->ymin + boxlist[i]->ymax);
    zcen = 0.5*(boxlist[i]->zmin + boxlist[i]->zmax);
    xcenterdata[i] = xcen;
    ycenterdata[i] = ycen;
    zcenterdata[i] = zcen;  
  }
  whichcut[0] = xcenterdata; whichcut[1] = ycenterdata; whichcut[2] = zcenterdata;
  
  //  Get the bounding box of the root node -- the entire set of bounding boxes.
  thisbox = getbox(0,npoly-1);
  icnt = 0;
  treebox[icnt] = thisbox;
  nextbranch[icnt] = icnt;
   
  max_tri_sequence_value = npoly-1;
  min_tri_sequence_value = 0;
 
  if ( max_tri_sequence_value == 0 ) return 1;  //  Just one triangle in the data set.

  //  Get the cutting direction for the next branch.  
  cuttingdir = getcuttingdirection(thisbox); 
     
  thistreestack = new TreeStack(min_tri_sequence_value,
	                                   max_tri_sequence_value,
					   cuttingdir,icnt);
  icnt++;                                         
  mytreestack.push_back(thistreestack);
    
  while ( mytreestack.size() > 0 ) {
    TreeStack* thisone = mytreestack[mytreestack.size()-1];
    mytreestack.pop_back();
    cuttingdir = thisone->cuttingdir;
    min_tri_sequence_value = thisone->min;
    max_tri_sequence_value = thisone->max;


    median_value = (min_tri_sequence_value + max_tri_sequence_value)/2;
    cutarray = whichcut[cuttingdir];
    
    find_the_median(median_value-min_tri_sequence_value,0,
                    max_tri_sequence_value-min_tri_sequence_value,
                    cutarray,&isequence[min_tri_sequence_value]);
    if ( min_tri_sequence_value == median_value ) {  //  This is a leaf.
//      FSBoundingBox *thisbox = boxlist[isequence[median_value]];
      thisbox = new FSBoundingBox(
                             boxlist[isequence[median_value]]->xmin,
			     boxlist[isequence[median_value]]->ymin, 
                             boxlist[isequence[median_value]]->zmin,
			     boxlist[isequence[median_value]]->xmax,
			     boxlist[isequence[median_value]]->ymax,
			     boxlist[isequence[median_value]]->zmax);


      treebox[icnt] = thisbox;
      nextbranch[icnt] = -isequence[median_value];
      nextbranch[thisone->sequence] = icnt;
      icnt++;
    } else {
      thisbox = getbox(min_tri_sequence_value,median_value);
        nextbranch[thisone->sequence] = icnt;
      treebox[icnt] = thisbox;
//      nextbranch[icnt] = 11111;
      cuttingdir = getcuttingdirection(thisbox);

      thistreestack = new TreeStack(min_tri_sequence_value,
	                                   median_value,cuttingdir,icnt);
      icnt++;                                    
      mytreestack.push_back(thistreestack);
    } 
    //  The delete that follows is because we need to delete items that
    //  were created and placed on mytreestack.  This one comes from the pop()
    
    delete thisone;
    
    if ( median_value+1 == max_tri_sequence_value ) {  //  This is a leaf.
//      FSBoundingBox *thisbox = boxlist[isequence[max_tri_sequence_value]];
      thisbox = new FSBoundingBox(
                             boxlist[isequence[max_tri_sequence_value]]->xmin,
			     boxlist[isequence[max_tri_sequence_value]]->ymin, 
                             boxlist[isequence[max_tri_sequence_value]]->zmin,
			     boxlist[isequence[max_tri_sequence_value]]->xmax,
			     boxlist[isequence[max_tri_sequence_value]]->ymax,
			     boxlist[isequence[max_tri_sequence_value]]->zmax);


      treebox[icnt] = thisbox;
      nextbranch[icnt] = -isequence[max_tri_sequence_value];
      icnt++;
    } else {
      thisbox = getbox(median_value+1,max_tri_sequence_value);
      treebox[icnt] = thisbox;
//      nextbranch[icnt] = 22222;
      cuttingdir = getcuttingdirection(thisbox);
      thistreestack = new TreeStack(median_value+1,
	                                   max_tri_sequence_value,
					   cuttingdir,icnt);
      icnt++;                                     
      mytreestack.push_back(thistreestack);    
    }     

  } 
 
  delete [] isequence;
  delete [] xcenterdata;
  delete [] ycenterdata;
  delete [] zcenterdata;
   
  return 1;
}

FSBoundingBox* KDTree::getbox(int min, int max)
{
FSBoundingBox *thisbox;
int i;
//  Makes a bounding box and sets its size to include all of the boxes in the
//  sequence of bounding boxes from min to max of index isequence[].  Returns pointer
//  to this box.

  thisbox = new FSBoundingBox(1.e20,1.e20,1.e20,-1.e20,-1.e20,-1.e20);
double xmin, ymin, zmin, xmax, ymax, zmax;
  
  for ( i = min; i <= max; i++ ) {
    xmin = boxlistptr[isequence[i]]->xmin;
    ymin = boxlistptr[isequence[i]]->ymin;
    zmin = boxlistptr[isequence[i]]->zmin;
    xmax = boxlistptr[isequence[i]]->xmax;
    ymax = boxlistptr[isequence[i]]->ymax;
    zmax = boxlistptr[isequence[i]]->zmax;

    if ( xmin < thisbox->xmin ) thisbox->xmin = xmin;
    if ( ymin < thisbox->ymin ) thisbox->ymin = ymin;
    if ( zmin < thisbox->zmin ) thisbox->zmin = zmin;    
    if ( xmax > thisbox->xmax ) thisbox->xmax = xmax;
    if ( ymax > thisbox->ymax ) thisbox->ymax = ymax;
    if ( zmax > thisbox->zmax ) thisbox->zmax = zmax;
     
  }

  return thisbox;
}

void KDTree::box_kdtree_intersect(FSBoundingBox& bbox, int *count, int *indexlist) const
{
int index, i, child;
std::stack<int> mystack;

  *count = 0;
  if ( (bbox.xmax < (treebox[0]->xmin - epsilonkd) ) ||
       (bbox.xmin > (treebox[0]->xmax + epsilonkd) ) ||
       (bbox.ymax < (treebox[0]->ymin - epsilonkd) ) ||
       (bbox.ymin > (treebox[0]->ymax + epsilonkd) ) ||
       (bbox.zmax < (treebox[0]->zmin - epsilonkd) ) ||
       (bbox.zmin > (treebox[0]->zmax + epsilonkd) ) )
    return;
  //  Gotta put something on the stack, so do the root node 
  //  evaluation here.
  if ( nextbranch[0] <= 0 ) {
    indexlist[*count] = -nextbranch[0]; *count += 1;
    return;
  }
  mystack.push(0);

  while ( mystack.size() > 0 ) {
    index = mystack.top();
    mystack.pop();
    child = nextbranch[index];
    
    for ( i = 0; i < 2; i++ ) {
      if ( (bbox.xmax < (treebox[child+i]->xmin - epsilonkd) ) ||
           (bbox.xmin > (treebox[child+i]->xmax + epsilonkd) ) ||
           (bbox.ymax < (treebox[child+i]->ymin - epsilonkd) ) ||
           (bbox.ymin > (treebox[child+i]->ymax + epsilonkd) ) ||
           (bbox.zmax < (treebox[child+i]->zmin - epsilonkd) ) ||
           (bbox.zmin > (treebox[child+i]->zmax + epsilonkd) ) )
      continue;
      if ( nextbranch[child+i] <= 0 ) {        
		indexlist[*count] = -nextbranch[child+i]; *count += 1;	
      } else {            
        mystack.push(child+i);   
      }  
    }
  }  
  return;
}


void KDTree::find_the_median(int k, int l, int r, double *array, int *ia)
{
int i, j;
double t;

  while ( r > l ) {
    t = array[ia[k]];
    i = l;
    j = r;

    SWAP(ia[l],ia[k]);
    if ( array[ia[r]] > t ) SWAP(ia[l],ia[r]);
    while ( i < j ) {
      SWAP(ia[i],ia[j]);
      i += 1; j -= 1;
      while ( array[ia[i]] < t ) i++;
      while ( array[ia[j]] > t ) j--;
    }
    if ( array[ia[l]] == t ) SWAP(ia[l],ia[j]);
    else {
      j += 1;
      SWAP(ia[r],ia[j]);
    }
    if ( j <= k ) l = j + 1;
    if ( k <= j ) r = j - 1;
  }

}

int KDTree::getcuttingdirection(FSBoundingBox* box)
{
double xlen, ylen, zlen;

  xlen = box->xmax - box->xmin;
  ylen = box->ymax - box->ymin;
  zlen = box->zmax - box->zmin;
  
  if ( (xlen >= ylen) && (xlen >= zlen) ) return XCUT; 
  else if ( ylen >= zlen ) return YCUT;
  else return ZCUT;
}


bool KDTree::rayintersectsbox(FSBoundingBox *box)
{
double pmin, pmax, tmin, tmax;
double dtemp;

//  Get the parametric distance along each direction from the min point to
//  the min and max box planes.  If the min dist is grater than the max
//  dist, we have to swap them.  Keep a running total of the max min and the
//  min max.  The ray intersects the box iff tmin <= tmax and tmax >= 0.0.
//  Otherwise, the ray misses the box or points away from the box (with the
//  starting point outside).

  tmin = (box->xmin - rayxstart)/dx;
  tmax = (box->xmax - rayxstart)/dx;
  if ( tmin > tmax ) {
  dtemp = tmin; tmin = tmax; tmax = dtemp;
  }
  pmin = (box->ymin - rayystart)/dy;
  pmax = (box->ymax - rayystart)/dy;
  if ( pmin > pmax ) {
  dtemp = pmin; pmin = pmax; pmax = dtemp;
  }
  tmin = MAXX(pmin,tmin);
  tmax = MINN(pmax,tmax);

  if ( tmin > tmax ) return false;
  
  pmin = (box->zmin - rayzstart)/dz;
  pmax = (box->zmax - rayzstart)/dz;
  if ( pmin > pmax ) {
  dtemp = pmin; pmin = pmax; pmax = dtemp;
  }
  tmin = MAXX(pmin,tmin);
  tmax = MINN(pmax,tmax);
  
  if ( (tmax < 0.0) || (tmin > tmax) ) return false;
  
  return true;
}
  
