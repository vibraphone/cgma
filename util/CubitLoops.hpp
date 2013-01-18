#ifndef CUBIT_LOOPS_HPP
#define CUBIT_LOOPS_HPP

#include "CubitDefines.h"
#include <vector>
#include <set>
#include <map>

template <class C, class V>
class CubitLoops
{
public:
  struct CoEdge
  {
    C* curve;
    V* start;
    V* end;
    CubitSense sense;
  };

  static bool make_loops(std::vector<CoEdge>& coedges,
      std::vector<std::vector<CoEdge*> >& loops);

private:

  static bool recursive_make_loop(V* start_vertex, CoEdge* coedge,
    std::set<CoEdge* >& used_coedges,
    std::multimap<V*, CoEdge*>& start_coedge_map,
    std::vector<CoEdge*>& loop);
};

template <class C, class V> inline
bool CubitLoops<C, V>::recursive_make_loop(V* start_vertex, CoEdge* coedge,
    std::set<CoEdge* >& used_coedges,
    std::multimap<V*, CoEdge*>& start_coedge_map,
    std::vector<CoEdge*>& loop)
{
  V* curr_vertex;
  if (coedge->sense == CUBIT_FORWARD)
    curr_vertex = coedge->end;
  else
    curr_vertex = coedge->start;
  loop.push_back(coedge);
  used_coedges.insert(coedge);

  while (curr_vertex != start_vertex) 
  {
    typename std::multimap<V*, CoEdge*>::iterator iter;
    typename std::multimap<V*, CoEdge*>::iterator last;

    iter = start_coedge_map.lower_bound(curr_vertex);
    last = start_coedge_map.upper_bound(curr_vertex);

    std::vector<CoEdge*> possible_coedges;
    for (/*preinitialized*/; iter != last; iter++)
    {
      if (used_coedges.find(iter->second) == used_coedges.end())
        possible_coedges.push_back(iter->second);
    }

    if (possible_coedges.size() == 0)
      return false;
    else if (possible_coedges.size() == 1)
    {
      coedge = possible_coedges[0];
      loop.push_back(coedge);
      used_coedges.insert(coedge);
      if (coedge->sense == CUBIT_FORWARD)
        curr_vertex = coedge->end;
      else
        curr_vertex = coedge->start;
    }
    else
    {
      for (size_t i=0; i<possible_coedges.size(); i++)
      {
        std::vector<CoEdge*> sub_loop;
        if (recursive_make_loop(curr_vertex, possible_coedges[i], used_coedges, start_coedge_map, sub_loop) )
        {
          loop.insert(loop.end(), sub_loop.begin(), sub_loop.end());
        }
        else
        {
          for (size_t j=0; j<sub_loop.size(); j++)
            used_coedges.erase(sub_loop[j]);
          coedge = possible_coedges[i];
        }
      }
      loop.push_back(coedge);
      used_coedges.insert(coedge);
      if (coedge->sense == CUBIT_FORWARD)
        curr_vertex = coedge->end;
      else
        curr_vertex = coedge->start;
    }
  }

  return true;
}

template <class C, class V> inline
bool CubitLoops<C, V>::make_loops(std::vector<CoEdge>& coedges,
      std::vector<std::vector<CoEdge*> >& loops)
{
  std::multimap<V*, CoEdge*> start_coedge_map;

  for (size_t i=0; i<coedges.size(); i++)
  {
    if (coedges[i].sense == CUBIT_FORWARD)
      start_coedge_map.insert(std::make_pair(coedges[i].start, &coedges[i]));
    else
      start_coedge_map.insert(std::make_pair(coedges[i].end, &coedges[i]));
  }

  std::set<CoEdge*> used_coedges;
  for (size_t i=0; i<coedges.size(); i++)
  {
    if (used_coedges.find(&coedges[i]) != used_coedges.end())
      continue;
    V* start_vertex;
    if (coedges[i].sense == CUBIT_FORWARD)
      start_vertex = coedges[i].start;
    else
      start_vertex = coedges[i].end;

    typename std::multimap<V*, CoEdge*>::iterator iter;
    typename std::multimap<V*, CoEdge*>::iterator last;
    iter = start_coedge_map.lower_bound(start_vertex);
    last = start_coedge_map.upper_bound(start_vertex);
    std::vector<CoEdge*> loop;

    for (/*preinitialized*/; iter != last; iter++)
    {
      std::vector<CoEdge*> sub_loop;
      recursive_make_loop(start_vertex, iter->second, used_coedges, 
                          start_coedge_map, loop);
    }

    if (loop.empty())
      return false;
    loops.push_back(std::vector<CoEdge*>(loop.begin(), loop.end()) );
  }

  if (coedges.size() != used_coedges.size())
    return false;

  return true;
}

#endif
