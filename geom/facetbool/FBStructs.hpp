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

#ifndef _FBSTRUCTS
#define _FBSTRUCTS
#include <vector>
#include <list>
#include "CubitDefines.h"
#include "FBDefines.hpp"

class FSBoundingBox {
public:
  FSBoundingBox() { }
  FSBoundingBox(float xmi, float ymi, float zmi, float xma, float yma, float zma)
  {
    xmin = xmi; xmax = xma; ymin = ymi; ymax = yma; zmin = zmi; zmax = zma;
  }
  ~FSBoundingBox() { }
  float  xmin, xmax, ymin, ymax, zmin, zmax;
};

typedef std::vector<FSBoundingBox* > FSBOXVECTOR;

class FB_Coord {

public:
   FB_Coord(double x, double y, double z) {
     coord[0] = x; coord[1] = y; coord[2] = z;
     is_on_boundary = false;
   }
   ~FB_Coord() {  }
   double coord[3];
   bool is_on_boundary;

private:

};

class FB_Edge {

public:
   FB_Edge(int vert0, int vert1, int vert0_type, int vert1_type, bool is_intersection) {
     v0 = vert0; v1 = vert1;
     v0_type = vert0_type; v1_type = vert1_type;
     mark = false;
     num_times = 0;
     slope = 0.0;
     quadrant = 0;
     if ( is_intersection == true ) 
       edge_type = INTERSECTION_EDGE;
     else 
       edge_type = UNKNOWN;
     is_intersection_edge = is_intersection;
     cubitedgeindex = 0;
   }
   FB_Edge();
   ~FB_Edge() {  }
   int v0, v1;
   int v0_type, v1_type;
   int edge_type;
   bool mark;
   bool is_intersection_edge;
   int num_times;
   double slope;       
   unsigned int quadrant;
   int cubitedgeindex;
  
      
  
  
      
};

class FB_Triangle {

public:
   FB_Triangle(int vert0, int vert1, int vert2, int myparent, int csvalue,
               int c0value, int c1value, int c2value) {
     v0 = vert0; v1 = vert1; v2 = vert2;
     dudded = false;
     parent = myparent;
     cubitsurfaceindex = csvalue;
     cubitedge0index = c0value;
     cubitedge1index = c1value;
     cubitedge2index = c2value;     
   }
   FB_Triangle();
   ~FB_Triangle() {
     std::vector<FB_Edge*>::iterator itt;
     itt = edge_list.begin();
     while ( itt != edge_list.end() ) {
       delete *itt;
       itt++;
     }
     edge_list.clear();
   }
    //! Util function to determine whether two verts are in a Triangle.
  bool are_fb_coords_in_fb_triangle(int fb_vert_index_1,
                                    int fb_vert_index_2)
    {
      if((fb_vert_index_1 == v0 || fb_vert_index_1 == v1 ||
          fb_vert_index_1 == v2) && (fb_vert_index_2== v0 ||
         fb_vert_index_2 == v1 || fb_vert_index_2 == v2)){
        return true;
      }
      return false;
    }

   int v0, v1, v2;
   FSBoundingBox boundingbox;
   double a, b, c, d; // plane coefficients
   bool dudded;
   std::vector<FB_Edge*> edge_list;
   int parent;
   int cubitsurfaceindex;
   int cubitedge0index, cubitedge1index, cubitedge2index;
};

class VertexStuff {

public:
  VertexStuff(int vert, int vert_type, double p_coord, double s_coord)
  {  v0 = vert;
     v0type = vert_type;
     p_dir_coord = p_coord;
     s_dir_coord = s_coord; 
  }
  
  ~VertexStuff() {
    edge_list.clear();
  }
  
  int v0;
  int v0type;
  double p_dir_coord, s_dir_coord;
  
  std::vector<FB_Edge*> edge_list;
};

class vertstuffcompfn_less {

public:
  bool operator()(const VertexStuff* vfirst, const VertexStuff* vsecond) const
  {
    if ( vfirst->p_dir_coord < vsecond->p_dir_coord ) return true;
    else if ( vfirst->p_dir_coord > vsecond->p_dir_coord ) return false;
    else {
      if ( vfirst->s_dir_coord < vsecond->s_dir_coord ) return true;
      else return false;
    }
  }
};

class edgecompfn_less {

public:
  bool operator()(const FB_Edge* efirst, const FB_Edge* esecond) const
  {
     if ( efirst->quadrant < esecond->quadrant ) return true;
     if ( efirst->quadrant > esecond->quadrant ) return false;

       if ( efirst->slope < esecond->slope ) return true;
       else return false;       

  }
};

class edgecompfn_more {

public:
  bool operator()(const FB_Edge* efirst, const FB_Edge* esecond) const
  {
     if ( efirst->quadrant > esecond->quadrant ) return true;
     if ( efirst->quadrant < esecond->quadrant ) return false;
//     if ( efirst->quadrant > esecond->quadrant ) return true;

       if ( efirst->slope < esecond->slope ) return false;
       else return true;       

  }
};

#endif
