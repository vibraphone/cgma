//- Class: ChollaVolume
//- Owner: Steven J. Owen
//- Description: volume representation for the Cholla entities
//- Created: 8/30/2009
//- Checked By:
//- Version:

#ifndef ChollaVolume_HPP
#define ChollaVolume_HPP

#include "DLIList.hpp"
#include "ChollaEntity.hpp"

class ChollaSurface;

class ChollaVolume : public ChollaEntity
{
private:   
  int id;
  int blockId;
  
  DLIList<ChollaSurface*> surfaceList;
  void *myVolume;
  
public:
   
  ChollaVolume(int block_id);
  ~ChollaVolume();

  void assign_geometric_volume(void *vol)
    {myVolume = vol;}
  void* get_geometric_volume()
    {return myVolume;}
  void get_surfaces( DLIList<ChollaSurface*> &cholla_surf_list )
    {cholla_surf_list = surfaceList; }
  void add_surface( ChollaSurface *cholla_surf_ptr )
    {surfaceList.append(cholla_surf_ptr);}
  void add_surface_unique( ChollaSurface *cholla_surf_ptr )
    {surfaceList.append_unique(cholla_surf_ptr);}
  void remove_surface( ChollaSurface *cholla_surf_ptr)
    {surfaceList.remove(cholla_surf_ptr);}
  int get_block_id()
    {return blockId;}
  void set_block_id(int flag)
    { blockId = flag; }

  int get_id(){return id;}

  void debug_draw();
};

#endif





