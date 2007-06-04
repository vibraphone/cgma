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
#include <algorithm>
#include "FBRetriangulate.hpp"
#include "FBTiler.hpp"
#include "CubitMessage.hpp"

#ifdef KEEP_BOYD10_KEEP
bool edgecf_less(const FB_Edge* efirst, const FB_Edge* esecond)
{
  if ( efirst->quadrant < esecond->quadrant ) return true;
  if ( efirst->quadrant > esecond->quadrant ) return false;

  if ( efirst->slope < esecond->slope ) return true;
  else return false;       

}
#endif
    
FBRetriangulate::FBRetriangulate(std::vector<FB_Coord *>& my_verts,
                                 std::vector<FB_Triangle *>& my_tris,
                                 std::vector<int>& my_newfacets,
                                 std::vector<int>& my_newfacetsindex)
{
  verts = my_verts;
  tris = &my_tris;
  newfacets = &my_newfacets;
  newfacetsindex = &my_newfacetsindex;
  p_dir=0;
  s_dir=0;
}
    
FBRetriangulate::FBRetriangulate(std::vector<FB_Coord *>& my_verts,
                                 std::vector<FB_Triangle *>& my_tris,
                                 std::vector<int>& my_newfacets)
{
  verts = my_verts;
  tris = &my_tris;
  newfacets = &my_newfacets;
  newfacetsindex = 0;
  p_dir=0;
  s_dir=0;
}

FBRetriangulate::~FBRetriangulate()
{

}

CubitStatus FBRetriangulate::retriangulate_this_tri(int sequence)
{
  CubitStatus status;
  double xspan, yspan, zspan;
  std::vector<FB_Triangle *>::iterator itt;

  itt = tris->begin();
  itt += sequence;
  my_tri = *itt;   
  status = CUBIT_SUCCESS;
  
//  tri = &itt;  
  xspan = my_tri->boundingbox.xmax - my_tri->boundingbox.xmin;
  yspan = my_tri->boundingbox.ymax - my_tri->boundingbox.ymin;
  zspan = my_tri->boundingbox.zmax - my_tri->boundingbox.zmin;
  if ( ( fabs(my_tri->c) >= fabs(my_tri->a) ) &&
       ( fabs(my_tri->c) >= fabs(my_tri->b) ) ) {
    p_dir = 1; 
    s_dir = 0;
    if ( xspan > yspan ) {
      p_dir = 0; 
      s_dir = 1;
    }
  }
  else if ( fabs(my_tri->b) >= fabs(my_tri->a) )  {
    p_dir = 2; 
    s_dir = 0;
    if ( xspan > zspan ) {
      p_dir = 0;
      s_dir = 2;
    }
  }
  else {
    p_dir = 1;
    s_dir = 2;
    if ( zspan > yspan ) {
      p_dir = 2;
      s_dir = 1;
    }
  }
   
  if ( p_dir == 0 ) {
    if ( s_dir == 1 ) {
      if ( my_tri->c > 0.0 )
        winding = CCW;
      else
        winding = CW;
    }
    else if (s_dir == 2 ) {
      if ( my_tri->b > 0.0 )
        winding = CW;
      else
        winding = CCW;
    }
  }
  else if ( p_dir == 1 ) {
    if ( s_dir == 0 ) {
      if ( my_tri->c > 0.0 )
        winding = CW;
      else
        winding = CCW;
    }
    else if (s_dir == 2 ) {
      if ( my_tri->a > 0.0 )
        winding = CW;
      else
        winding = CCW;
    }
  }
  else if ( p_dir == 2 ) {
    if ( s_dir == 0 ) {
      if ( my_tri->b > 0.0 )
        winding = CCW;
      else
        winding = CW;
    }
    else if (s_dir == 1 ) {
      if ( my_tri->a > 0.0 )
        winding = CCW;
      else
        winding = CW;
    } 
  }
  else{
    PRINT_ERROR("Unexpected result.\n");
    return CUBIT_FAILURE;
  }
  
//if ( winding == CCW ) winding = CW;
//else winding = CCW;
  classify_edges();
  
//  Add edges around the perimeter of the triangle to its edge list.
  status = add_bdry_edges();
  make_vert_list();
  
//  Add internal edges where there is a local min or max 
  status = remove_min_max();

// Get the vertex chains that will be retriangulated.  
  
  bool chain_status = true; 
  sort_vertstufflist_edges();  
  std::vector<int> *chainlist = new std::vector<int>;
//  std::vector<FB_Triangle* > new_tris;
  unsigned int i, nfsize;
  nfsize = newfacets->size();
  while (chain_status == true ) {
    chain_status = get_a_chain(chainlist);
    if ( chain_status == false ) break;
    FBTiler *tiler = new FBTiler(verts,p_dir,s_dir,sequence,
                                 my_tri->a,my_tri->b,my_tri->c,newfacets);
    status = tiler->Tile_Poly(chainlist);
    delete tiler;
    chainlist->clear();

  }

  delete chainlist;
  unsigned int number_new;
  int e0index, e1index, e2index;

  number_new = newfacets->size() - nfsize;
  if ( number_new > 0 ) {
    std::vector<int>::iterator itp;
    itp = newfacets->begin();
    itp += nfsize;
    for ( i = 0; i < number_new; i += 3 ) {
      FB_Triangle *new_tri;
      e0index = e1index = e2index = 0;
      get_edge_indices(*itp,*(itp+1),*(itp+2),sequence,
                       e0index, e1index, e2index);
      new_tri = new FB_Triangle(*itp,*(itp+1),*(itp+2),
                                sequence,my_tri->cubitsurfaceindex,
                                e0index, e1index, e2index);
      tris->push_back(new_tri);                         
      itp += 3;   
    }
    if ( newfacetsindex ) newfacetsindex->push_back(nfsize);
  }  
  for ( i = 0; i < vertstufflist.size(); i++ ) {
    delete vertstufflist[i];
  }

  vertstufflist.clear();

  return status;

}

CubitStatus FBRetriangulate::remove_min_max()
{
  CubitStatus status;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe;
  bool goes_down, goes_up;
  VertexStuff *vstuff;
  unsigned int i;
  int this_vert, that_vert;
  double this_vert_p_dir_coord, that_vert_p_dir_coord;
  
  status = CUBIT_SUCCESS;

  for ( i = 0; i < vertstufflist.size(); i++ ) {
    vstuff = vertstufflist[i];
    if ( vstuff->v0type != INTERIOR_VERT ) continue;
    goes_down = goes_up = false;
    this_vert = vstuff->v0;
    this_vert_p_dir_coord = verts[this_vert]->coord[p_dir];
    dpe = vstuff->edge_list.begin();
    while ( dpe != vstuff->edge_list.end() ) {
      edge = *dpe;
      dpe++;
      if ( edge->v0 == this_vert ) that_vert = edge->v1;
      else that_vert = edge->v0;
      that_vert_p_dir_coord = verts[that_vert]->coord[p_dir]; 
      if ( that_vert_p_dir_coord < this_vert_p_dir_coord )
        goes_down = true;
      else if ( that_vert_p_dir_coord > this_vert_p_dir_coord )
        goes_up = true;
      else {
        double this_vert_s_dir_coord, that_vert_s_dir_coord;
        this_vert_s_dir_coord = verts[this_vert]->coord[s_dir];
        that_vert_s_dir_coord = verts[that_vert]->coord[s_dir]; 
        if ( this_vert_s_dir_coord < that_vert_s_dir_coord )
          goes_up = true;
        else goes_down = true;
      }
    }

    if ( goes_down == false ) {
      add_edge_down(this_vert,i);
    } else if ( goes_up == false ) {
      add_edge_up(this_vert,i);    
    } 
  }
         
  return status;
  
}

CubitStatus FBRetriangulate::add_bdry_edges()
{
//  For each triangle edge, make a list of pairs of edge points and
//  vertex numbers that touch the edge.  Then sort the list by distance
//  from the start of the edge.  (This was what we really put in the
//  pair, not the point itself.)  Finally, go through this sorted list
//  and add edges to the edge_list.  If the edge_list is empty, we
//  still must add in the triangle edge itself.

  CubitStatus status;
  FB_Edge *edge;
  std::list< std::pair<double,int> > edge0_list, edge1_list, edge2_list;
  std::list< std::pair<double,int> >::iterator dp;
  std::vector<FB_Edge*>::iterator dpe, dpe_orig;
  double vx0, vy0, vz0, vx1, vy1, vz1, vx2, vy2, vz2, dist;
  std::pair<double,int> mypair;

  status = CUBIT_SUCCESS;

  vx0 = verts[my_tri->v0]->coord[0];
  vy0 = verts[my_tri->v0]->coord[1];
  vz0 = verts[my_tri->v0]->coord[2];
  vx1 = verts[my_tri->v1]->coord[0];
  vy1 = verts[my_tri->v1]->coord[1];
  vz1 = verts[my_tri->v1]->coord[2];
  vx2 = verts[my_tri->v2]->coord[0];
  vy2 = verts[my_tri->v2]->coord[1];
  vz2 = verts[my_tri->v2]->coord[2];

  dpe_orig = my_tri->edge_list.end();
  dpe = my_tri->edge_list.begin();
  while ( dpe != my_tri->edge_list.end() ) {
    edge = *dpe;
    dpe++;
    if ( (edge->v0_type == INTERIOR_VERT) &&
         (edge->v1_type == INTERIOR_VERT) )
      continue;
    if ( edge->v0_type == EDGE_0 ) {
      dist = get_dist(vx0,vy0,vz0,verts[edge->v0]->coord[0],
                      verts[edge->v0]->coord[1],verts[edge->v0]->coord[2]);
      edge0_list.push_back(std::pair<double,int>(dist,edge->v0));
    } else if ( edge->v0_type == EDGE_1 ) { 
      dist = get_dist(vx1,vy1,vz1,verts[edge->v0]->coord[0],
                      verts[edge->v0]->coord[1],verts[edge->v0]->coord[2]);
      edge1_list.push_back(std::pair<double,int>(dist,edge->v0));
    } else if ( edge->v0_type == EDGE_2 ) {
      dist = get_dist(vx2,vy2,vz2,verts[edge->v0]->coord[0],
                      verts[edge->v0]->coord[1],verts[edge->v0]->coord[2]);
      edge2_list.push_back(std::pair<double,int>(dist,edge->v0));
    } 
    if ( edge->v1_type == EDGE_0 ) {
      dist = get_dist(vx0,vy0,vz0,verts[edge->v1]->coord[0],
                      verts[edge->v1]->coord[1],verts[edge->v1]->coord[2]);
      edge0_list.push_back(std::pair<double,int>(dist,edge->v1));
    } else if ( edge->v1_type == EDGE_1 ) {
      dist = get_dist(vx1,vy1,vz1,verts[edge->v1]->coord[0],
                      verts[edge->v1]->coord[1],verts[edge->v1]->coord[2]);
      edge1_list.push_back(std::pair<double,int>(dist,edge->v1));
    } else if ( edge->v1_type == EDGE_2 ) {
      dist = get_dist(vx2,vy2,vz2,verts[edge->v1]->coord[0],
                      verts[edge->v1]->coord[1],verts[edge->v1]->coord[2]);
      edge2_list.push_back(std::pair<double,int>(dist,edge->v1));
    }     
  }
  
  edge0_list.sort(); 
  edge1_list.sort(); 
  edge2_list.sort();
   
//  Now we have to remove all BDRY_EDGEs because they will be made anew
//  in what follows.  Erasing elements from a vector is inefficient, but
//  this shouldn't happen often.
  dpe = my_tri->edge_list.begin();

  while ( dpe != dpe_orig ) {
    edge = *dpe;
    if ( edge->edge_type == BDRY_EDGE ) {
      my_tri->edge_list.erase(dpe);
      dpe_orig -= 1;
    } else {
      dpe++;
    }
  }  
  int newv0, newv1, newv0_type, newv1_type;
  newv0_type = newv1_type = EDGE_0;
  newv0 = my_tri->v0;  
  dp = edge0_list.begin();
  while ( dp != edge0_list.end() ) {
    mypair = *dp;
    newv1 = mypair.second;
      //  It is possible to get the same vert more than once in the list.
      //  After sorting, they will be adjacent.  The following if statement
      //  causes duplicate verts to be used only once.
    if ( newv0 != newv1 ) {
      add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type);
      newv0 = newv1;
    }
    dp++;
  }
  newv1 = my_tri->v1;
  if ( newv0 != newv1 )
    add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type);

  newv0_type = newv1_type = EDGE_1;
  newv0 = my_tri->v1;  
  dp = edge1_list.begin();
  while ( dp != edge1_list.end() ) {
    mypair = *dp;
    newv1 = mypair.second;
    if ( newv0 != newv1 ) {
      add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type);
      newv0 = newv1;
    }
    dp++;
  }  
  newv1 = my_tri->v2;
  if ( newv0 != newv1 )
    add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type);

  newv0_type = newv1_type = EDGE_2;
  newv0 = my_tri->v2;  
  dp = edge2_list.begin();
  while ( dp != edge2_list.end() ) {
    mypair = *dp;
    newv1 = mypair.second;
    if ( newv0 != newv1 ) {
      add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type);
      newv0 = newv1;
    }
    dp++;
  }
  newv1 = my_tri->v0;
  if ( newv0 != newv1 )
    add_this_bdry_edge(newv0,newv1,newv0_type,newv1_type); 
    
  return status;
  
}

void FBRetriangulate::add_this_bdry_edge(int v0, int v1, int v0_type,
                                         int v1_type) 
{
//  Add an edge if it doesn't already exist.
  std::vector<FB_Edge*>::iterator dpe;
  bool ifoundit;
  FB_Edge *edge, *new_edge;

  ifoundit = false;
  dpe = my_tri->edge_list.begin();
  while ( dpe != my_tri->edge_list.end() ) {
    edge = *dpe;
    if ( ( ((int)edge->v0 == v0) && ((int)edge->v1 == v1) ) ||
         ( ((int)edge->v0 == v1) && ((int)edge->v1 == v0) ) ) {
      ifoundit = true; 
      break;
    }
    dpe++;
  }
  if ( ifoundit == false ) {
    new_edge = new FB_Edge(v0,v1,v0_type,v1_type,true);
    new_edge->edge_type = BDRY_EDGE;
    my_tri->edge_list.push_back(new_edge);  
  }    
}

void FBRetriangulate::make_vert_list()
{
  unsigned int i;
  int v0, v1;
  FB_Edge *edge;
  VertexStuff *vstuff;
  std::vector<FB_Edge*>::iterator dpe;
  bool foundv0, foundv1; 

  dpe = my_tri->edge_list.begin();
  while ( dpe != my_tri->edge_list.end() ) {
    edge = *dpe;
    dpe++;
    v0 = edge->v0;
    v1 = edge->v1;
    foundv0 = foundv1 = false;
    
    for ( i = 0; i < vertstufflist.size(); i++ ) {
      if ( vertstufflist[i]->v0 == v0 ) {
        vertstufflist[i]->edge_list.push_back(edge);
        foundv0 = true;
      }
      if ( vertstufflist[i]->v0 == v1 ) {
        vertstufflist[i]->edge_list.push_back(edge);
        foundv1 = true;
      }
    }
    if ( foundv0 == false ) {
      vstuff = new VertexStuff(v0, edge->v0_type,verts[v0]->coord[p_dir],
                               verts[v0]->coord[s_dir]);
      vstuff->edge_list.push_back(edge);
      vertstufflist.push_back(vstuff);
    }
    if ( foundv1 == false ) {
      vstuff = new VertexStuff(v1, edge->v1_type,verts[v1]->coord[p_dir],
                               verts[v1]->coord[s_dir]);
      vstuff->edge_list.push_back(edge);
      vertstufflist.push_back(vstuff);
    }
    
  } 
  
  std::vector<VertexStuff* >::iterator vitbegin, vitend;
  vitbegin = vertstufflist.begin();
  vitend = vertstufflist.end();

  std::sort(vitbegin,vitend,vertstuffcompfn_less());  

}

void FBRetriangulate::add_edge_up(int v0, int seq)
{
  unsigned int i, k;
  int v1, v0test, v1test, v0type, v1type;
  bool foundit, crossed;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe;

  v1 = 0; //To keep compiler from warning that v1 might be used uninitialized.
  v1type = 0;
  v0type = vertstufflist[seq]->v0type;
  for ( k = seq+1; k < vertstufflist.size(); k++ ) {
    foundit = true;
    v1 = vertstufflist[k]->v0;
    if ( fabs(verts[v1]->coord[p_dir]-verts[v0]->coord[p_dir]) < EPSILON )
      continue;
    v1type = vertstufflist[k]->v0type;
      //  v0 to v1 is the putative new edge.  Test whether it crosses any
      //  existing internal edge.  Any such internal edge that it crosses
      //  has to have an endpoint higher than v1.  If v1 is the next-to-top
      //  vertex, use it.
    if ( k == vertstufflist.size() - 1 ) {
      foundit = true;
      add_tri_edge(v0,v1,v0type,v1type);
      return;      
    }
    for ( i = k+1; i < vertstufflist.size()-1; i++ ) {
      dpe = vertstufflist[i]->edge_list.begin();

      while ( dpe != vertstufflist[i]->edge_list.end() ) {
        edge = *dpe;
        dpe++;
        if ( (edge->v0_type != INTERIOR_VERT) &&
             (edge->v1_type != INTERIOR_VERT) )
          continue;
        v0test = edge->v0;
        v1test = edge->v1;
        if ( (v0test == v1) || (v1test == v1) ) continue;
        crossed = test_for_crossing(v0,v1,v0test,v1test);
        if ( crossed == true ) {
          foundit = false;
          break;
        }
      }
      if ( foundit == false ) break;
    }  //  end of "for (i ...."
    if ( foundit == true ) break;
  }
  add_tri_edge(v0,v1,v0type,v1type);
}

void FBRetriangulate::add_edge_down(int v0, int seq)
{
  int i, k;
  int v1, v0test, v1test, v0type, v1type;
  bool foundit, crossed;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe;

  v1 = 0;  //To keep compiler from warning that v1 might be used uninitialized.
  v1type = 0;
  v0type = vertstufflist[seq]->v0type;
  for ( k = seq-1; k > -1; k-- ) {
    foundit = true;
    v1 = vertstufflist[k]->v0;
    if ( fabs(verts[v1]->coord[p_dir]-verts[v0]->coord[p_dir]) < EPSILON ) continue;
    v1type = vertstufflist[k]->v0type;
      //  v0 to v1 is the putative new edge.  Test whether it crosses
      //  any existing internal edge.  Any such internal edge that it
      //  crosses has to have an endpoint higher than v1.  If v1 is the
      //  next-to-top vertex, use it.
    if ( k == 0 ) {
      foundit = true;
      add_tri_edge(v0,v1,v0type,v1type);
      return;      
    }    
    for ( i = k-1; i > 0; i-- ) {
      dpe = vertstufflist[i]->edge_list.begin();

      while ( dpe != vertstufflist[i]->edge_list.end() ) {
        edge = *dpe;
        dpe++;
        if ( (edge->v0_type != INTERIOR_VERT) &&
             (edge->v1_type != INTERIOR_VERT) )
          continue;
        v0test = edge->v0;
        v1test = edge->v1;
        if ( (v0test == v1) || (v1test == v1) ) continue;
        crossed = test_for_crossing(v0,v1,v0test,v1test);
        if ( crossed == true ) {
          foundit = false;
          break;
        }
      }
      if ( foundit == false ) break;
    }  //  end of "for (i ...."
    if ( foundit == true ) break;
  }
  add_tri_edge(v0,v1,v0type,v1type);
}
  
bool FBRetriangulate::test_for_crossing(int v0, int v1, int v2, int v3)
{
  double x0, y0, x1, y1, x2, y2, x3, y3, dxa, dya, dxb, dyb, p01x, p01y;
  double product, dasq, dbsq, prodsq;
  double s, t;

  x0 = verts[v0]->coord[p_dir]; y0 = verts[v0]->coord[s_dir];
  x1 = verts[v1]->coord[p_dir]; y1 = verts[v1]->coord[s_dir];
  x2 = verts[v2]->coord[p_dir]; y2 = verts[v2]->coord[s_dir];
  x3 = verts[v3]->coord[p_dir]; y3 = verts[v3]->coord[s_dir];
  dxa = x1 - x0; dya = y1 - y0;
  dxb = x3 - x2; dyb = y3 - y2;
  
  product = dxa*dyb - dya*dxb;
  dasq = dxa*dxa + dya*dya;
  dbsq = dxb*dxb + dyb*dyb;
  prodsq = product*product;
  
  if ( prodsq > EPSILON2*dasq*dbsq ) {
    p01x = x2 - x0;
    p01y = y2 - y0;
    s = (p01x*dyb - p01y*dxb)/product;
    if ( (s < 0.0) || (s > 1.0) ) return false; 
  
    t = (p01x*dya - p01y*dxa)/product;
    if ( (t < 0.0) || (t > 1.0) ) return false; 
    
  }
  
  return true;
}

void FBRetriangulate::add_tri_edge(int v0, int v1, int v0_type, int v1_type)
{
  FB_Edge *edge;
  unsigned int i;

  edge = new FB_Edge(v0,v1,v0_type,v1_type,false);
  edge->edge_type = INTERIOR_EDGE;
  my_tri->edge_list.push_back(edge);
  for ( i = 0; i < vertstufflist.size(); i++ ) {
    if ( vertstufflist[i]->v0 == v0 ) {
      vertstufflist[i]->edge_list.push_back(edge);
    }
    if ( vertstufflist[i]->v0 == v1 ) {
      vertstufflist[i]->edge_list.push_back(edge);
      v1_type = vertstufflist[i]->v0type;
    }
  }
  
}

bool FBRetriangulate::get_a_chain(std::vector<int> *chainlist)
{
  bool status;
  int vthis, vprev, vstart;
  int direction;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe;
  status = false;
  dpe = my_tri->edge_list.begin();
    // to keep compiler from warning that edge might be used uninitialized
  edge = 0;
  
  while ( dpe !=  my_tri->edge_list.end() ) {
    edge = *dpe;
    dpe++;
      //  Skip edges that are not interior edges.
//    if ( edge->v0_type == edge->v1_type ) continue;
    if ( edge->edge_type == BDRY_EDGE  ) continue;
    if ( edge->num_times < 2 ) {
      if ( (edge->edge_type == BDRY_EDGE) && (edge->num_times == 1) )
        continue;
      status = true;
      edge->num_times++;
      break;
    }
  }
  if ( status == true ) {
    if ( edge->num_times == 1 ) {
      vthis = edge->v1;
      vprev = edge->v0;
      direction = 1;
    } else {
      vthis = edge->v0;
      vprev = edge->v1;
      direction = 1;
    }
    vstart = vprev;
    chainlist->push_back(vthis);
    
    while ( vthis != vstart ) {
      get_next_vert_in_chain(vthis,vprev,direction);
      chainlist->push_back(vthis);
    }    
  }

  if ( status == false ) { // There were no interior edges.
                           // So just go around the perimeter.
    dpe = my_tri->edge_list.begin();
    edge = *dpe;
    if ( edge->num_times != 0 ) {
      status = false;
      return status;
    }
    edge->num_times = 1;
    vthis = edge->v0;
    vprev = edge->v1;
    direction = 1;
    vstart = vprev;
    chainlist->push_back(vthis);
    
    while ( vthis != vstart ) {
      get_next_vert_in_chain(vthis,vprev,direction);
      chainlist->push_back(vthis);
    } 
    status = true;  
  }
  
  return status;
 
}

void FBRetriangulate::get_next_vert_in_chain(int& vthis,
                                             int& vprev, int direction)
//  Right now, direction always equals 1.
{
  unsigned int i;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe, dpebegin, dpeend;

  for ( i = 0; i < vertstufflist.size(); i++ ) {
    if ( vertstufflist[i]->v0 == (int)vthis ) break;
  }
  dpe = dpebegin = vertstufflist[i]->edge_list.begin();
  dpeend = vertstufflist[i]->edge_list.end();
  while ( dpe != dpeend ) {
    edge = *dpe;
    if ( (edge->v0 == vprev) || (edge->v1 == vprev) ) break;  
    dpe++;
  }
  if ( direction == 1 ) {
    if ( dpe == dpebegin )
      dpe = vertstufflist[i]->edge_list.end();
    dpe--;    
  } else {
    dpe++;
    if ( dpe == dpeend )
      dpe = dpebegin;
  }
  edge = *dpe;
  vprev = vthis;
  if ( vthis == edge->v1 ) {
      //  Swap the edge verts.  We do this so that next time we see this
      //  edge we will know that it was oriented to point in the direction
      //  of the previous loop.  Thus we will know to proceed in the
      //  opposite direction next time.
      //  (see edge->num_times test in get_a_chain().)
    unsigned int vtemp;
    vtemp = edge->v0;
    edge->v0 = edge->v1;
    edge->v1 = vtemp;
    vtemp = edge->v0_type;
    edge->v0_type = edge->v1_type;
    edge->v1_type = vtemp;
  }
  vthis = edge->v1;
  edge->num_times++;
  
}

void FBRetriangulate::sort_vertstufflist_edges()
{
//  For each vertex that has more than 2 edges, sort the edges in CCW order.
  unsigned int i;
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe, dpebegin, dpeend;
  double x0, y0, x1, y1, slope;
  int v0;
  unsigned int quadrant;

  quadrant = 0;  //  Initialize so compiler won't warn.
  
  for ( i = 0; i < vertstufflist.size(); i++ ) {
    if ( vertstufflist[i]->edge_list.size() > 2 ) {
      v0 = vertstufflist[i]->v0;
      x0 = verts[v0]->coord[p_dir];
      y0 = verts[v0]->coord[s_dir];      
      dpe = vertstufflist[i]->edge_list.begin();
      while ( dpe != vertstufflist[i]->edge_list.end() ) {
        edge = *dpe;
        dpe++;
        if ( edge->v0 == v0 ) {
          x1 = verts[edge->v1]->coord[p_dir];
          y1 = verts[edge->v1]->coord[s_dir];
        } else {
          x1 = verts[edge->v0]->coord[p_dir];
          y1 = verts[edge->v0]->coord[s_dir];
        }
        if ( fabs(x1-x0) < EPSILON ) {
          if ( y1 > y0 ) {
            slope = CUBIT_DBL_MAX;
            quadrant = 1;
          } else {
            slope = -CUBIT_DBL_MAX;
            quadrant = 4;
          }
        } else {
          slope = (y1-y0)/(x1-x0);
          if ( (x1 >= x0) && (y1 >= y0) ) quadrant = 1;
          else if ( (x1 < x0) && (y1 > y0) ) quadrant = 2;
          else if ( (x1 <= x0) && (y1 <= y0) ) quadrant = 3;
          else if ( (x1 > x0) && (y1 < y0) ) quadrant = 4;
        }
        edge->slope = slope;
        edge->quadrant = quadrant;
      }
        //  Now sort the edge list by the value of quadrant or slope.
//      vertstufflist[i]->edge_list.sort(edgecompfn_less());
//      vertstufflist[i]->edge_list.sort(edgecf_less);
      dpebegin = vertstufflist[i]->edge_list.begin();
      dpeend = vertstufflist[i]->edge_list.end();   
//      if ( winding == CW )   
      std::sort(dpebegin,dpeend,edgecompfn_less());
//      else
//        std::sort(dpebegin,dpeend,edgecompfn_more());
      
    }

    dpe = vertstufflist[i]->edge_list.begin();
    while ( dpe != vertstufflist[i]->edge_list.end() ) {
      edge = *dpe;
      dpe++;
    }
  
  
  }
}

void FBRetriangulate::classify_edges()
{
//  Flag boundary edges as such.
  FB_Edge *edge;
  std::vector<FB_Edge*>::iterator dpe;
  int type;

  dpe = my_tri->edge_list.begin();
  while ( dpe != my_tri->edge_list.end() ) {
    edge = *dpe;
    
    if ( (edge->v0_type != INTERIOR_VERT) &&
         (edge->v1_type != INTERIOR_VERT) ) {
      if ( edge->v0_type == edge->v1_type )
        edge->edge_type = BDRY_EDGE;
      else {
        type = UNKNOWN;
        switch( edge->v0_type ) {
          case VERTEX_0:
            if ( edge->v1_type != EDGE_1 ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE) {
                if ( edge->v1_type == VERTEX_1 ) 
                  my_tri->cubitedge0index = INTERSECTION_EDGE;
                else if ( edge->v1_type == VERTEX_2 ) 
                  my_tri->cubitedge2index = INTERSECTION_EDGE;
              }
            }          
            break;
          case VERTEX_1:
            if ( edge->v1_type != EDGE_2 ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE) {
                if ( edge->v1_type == VERTEX_2 ) 
                  my_tri->cubitedge1index = INTERSECTION_EDGE;
                else if ( edge->v1_type == VERTEX_0 ) 
                  my_tri->cubitedge0index = INTERSECTION_EDGE;
              } 
            }          
            break;
          case VERTEX_2:
            if ( edge->v1_type != EDGE_0 ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE) {
                if ( edge->v1_type == VERTEX_1 ) 
                  my_tri->cubitedge1index = INTERSECTION_EDGE;
                else if ( edge->v1_type == VERTEX_0 ) 
                  my_tri->cubitedge2index = INTERSECTION_EDGE;
              }
            }          
            break;
          case EDGE_0:
            if ( (edge->v1_type == VERTEX_0) ||
                 (edge->v1_type == VERTEX_1) ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE)
                my_tri->cubitedge0index = INTERSECTION_EDGE;
            }
            break;
          case EDGE_1:
            if ( (edge->v1_type == VERTEX_1) ||
                 (edge->v1_type == VERTEX_2) ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE)
                my_tri->cubitedge1index = INTERSECTION_EDGE;
            }
            break;
          case EDGE_2:
            if ( (edge->v1_type == VERTEX_2) ||
                 (edge->v1_type == VERTEX_0) ) {
              type = BDRY_EDGE;
              if ( edge->edge_type == INTERSECTION_EDGE)
                my_tri->cubitedge2index = INTERSECTION_EDGE;
            }
            break;
        }
        if ( type == BDRY_EDGE ) edge->edge_type = BDRY_EDGE;      
      }        
    } // else edge->edge_type = INTERIOR_EDGE;
    dpe++;
  }
       
}

void FBRetriangulate::get_edge_indices(int v0, int v1, int v2, int parent, 
                                       int &e0index, int &e1index,
                                       int &e2index)
{
  std::vector<FB_Edge*>::iterator itt;
  std::vector<FB_Triangle*>::iterator itp;
  FB_Edge *edge;
  int e_v0, e_v1, e_type, e_v0type, e_v1type;

  itp = tris->begin();  
  itp += parent;
  itt = (*itp)->edge_list.begin();
  while ( itt != (*itp)->edge_list.end() ) {
    edge = *itt;
    e_v0 = edge->v0;
    e_v1 = edge->v1;
    e_v0type = edge->v0_type;
    e_v1type = edge->v1_type;
    e_type = edge->edge_type;
    if ( ( (e_v0 == v0) && (e_v1 == v1) ) ||
         ( (e_v0 == v1) && (e_v1 == v0) ) ) {
      if ( e_type == INTERSECTION_EDGE )
        e0index = INTERSECTION_EDGE; 
      else if ( e_type == INTERIOR_EDGE )
        e0index = 0;  
      else if ( (e_v0type == EDGE_0) || (e_v1type == EDGE_0) )
        e0index = (*itp)->cubitedge0index;
      else if ( (e_v0type == EDGE_1) || (e_v1type == EDGE_1) )
        e0index = (*itp)->cubitedge1index;
      else if ( (e_v0type == EDGE_2) || (e_v1type == EDGE_2) )
        e0index = (*itp)->cubitedge2index;      
      else if ( ( (e_v0type == VERTEX_0) || (e_v0type == VERTEX_1) ) &
                ( (e_v1type == VERTEX_0) || (e_v1type == VERTEX_1) ) )
        e0index = (*itp)->cubitedge0index;
      else if ( ( (e_v0type == VERTEX_1) || (e_v0type == VERTEX_2) ) &
                ( (e_v1type == VERTEX_1) || (e_v1type == VERTEX_2) ) )
        e1index = (*itp)->cubitedge1index;
      else if ( ( (e_v0type == VERTEX_2) || (e_v0type == VERTEX_0) ) &
                ( (e_v1type == VERTEX_2) || (e_v1type == VERTEX_0) ) )
        e2index = (*itp)->cubitedge2index;                
    }
    else if ( ( (e_v0 == v1) && (e_v1 == v2) ) ||
              ( (e_v0 == v2) && (e_v1 == v1) ) ) {
      if ( e_type == INTERSECTION_EDGE )
        e1index = INTERSECTION_EDGE; 
      else if ( e_type == INTERIOR_EDGE )
        e1index = 0;    
      else if ( (e_v0type == EDGE_0) || (e_v1type == EDGE_0) )
        e1index = (*itp)->cubitedge0index;
      else if ( (e_v0type == EDGE_1) || (e_v1type == EDGE_1) )
        e1index = (*itp)->cubitedge1index;
      else if ( (e_v0type == EDGE_2) || (e_v1type == EDGE_2) )
        e1index = (*itp)->cubitedge2index;
      else if ( ( (e_v0type == VERTEX_0) || (e_v0type == VERTEX_1) ) &
                ( (e_v1type == VERTEX_0) || (e_v1type == VERTEX_1) ) )
        e0index = (*itp)->cubitedge0index;
      else if ( ( (e_v0type == VERTEX_1) || (e_v0type == VERTEX_2) ) &
                ( (e_v1type == VERTEX_1) || (e_v1type == VERTEX_2) ) )
        e1index = (*itp)->cubitedge1index;
      else if ( ( (e_v0type == VERTEX_2) || (e_v0type == VERTEX_0) ) &
                ( (e_v1type == VERTEX_2) || (e_v1type == VERTEX_0) ) )
        e2index = (*itp)->cubitedge2index;
    }
    else if ( ( (e_v0 == v2) && (e_v1 == v0) ) ||
              ( (e_v0 == v0) && (e_v1 == v2) ) ) {
      if ( e_type == INTERSECTION_EDGE )
        e2index = INTERSECTION_EDGE; 
      else if ( e_type == INTERIOR_EDGE )
        e2index = 0;    
      else if ( (e_v0type == EDGE_0) || (e_v1type == EDGE_0) )
        e2index = (*itp)->cubitedge0index;
      else if ( (e_v0type == EDGE_1) || (e_v1type == EDGE_1) )
        e2index = (*itp)->cubitedge1index;
      else if ( (e_v0type == EDGE_2) || (e_v1type == EDGE_2) )
        e2index = (*itp)->cubitedge2index;
      else if ( ( (e_v0type == VERTEX_0) || (e_v0type == VERTEX_1) ) &
                ( (e_v1type == VERTEX_0) || (e_v1type == VERTEX_1) ) )
        e0index = (*itp)->cubitedge0index;
      else if ( ( (e_v0type == VERTEX_1) || (e_v0type == VERTEX_2) ) &
                ( (e_v1type == VERTEX_1) || (e_v1type == VERTEX_2) ) )
        e1index = (*itp)->cubitedge1index;
      else if ( ( (e_v0type == VERTEX_2) || (e_v0type == VERTEX_0) ) &
                ( (e_v1type == VERTEX_2) || (e_v1type == VERTEX_0) ) )
        e2index = (*itp)->cubitedge2index;
    }
    
    itt++;
  }

}

                                       
