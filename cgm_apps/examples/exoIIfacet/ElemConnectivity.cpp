#include "ElemConnectivity.hpp"
#include "exodusII.h"
#include <string.h>
#include "CubitMessage.hpp"

const char* const ElemConnectivity::element_type_names[] = {
  "BAR",   "BAR2",   "BAR3", 
  "BEAM",  "BEAM2",  "BEAM3",
  "TRUSS", "TRUSS2", "TRUSS3",

  "TRI",   "TRI3",   "TRI4",   "TRI6",   "TRI7",
  "SHELL", "SHELL4", "SHELL5", "SHELL8", "SHELL9",
  "QUAD",  "QUAD4",  "QUAD5",  "QUAD8",  "QUAD9",

  "TETRA",  "TETRA4",  "TETRA5",  "TETRA8",   "TETRA10",  "TETRA14",  "TETRA15",
  "PYRAMID","PYRAMID5","PYRAMID6","PYRAMID10","PYRAMID13","PYRAMID18","PYRAMID19",
  "WEDGE",  "WEDGE6",  "WEDGE7",  "WEDGE11",  "WEDGE15",  "WEDGE20",  "WEDGE21",
  "HEX",    "HEX8",    "HEX9",    "HEX14",    "HEX20",    "HEX26",    "HEX27", 

  "HEXSHELL","UNKNOWN_ELEMENT_TYPE" };


ElemConnectivity::ElemConnectivity( int exo_file_id, int num_exo_blocks )
  : file_id     ( exo_file_id       ), 
    block_count ( num_exo_blocks    ),
    max_block_id( num_exo_blocks + 1),
    last_block  ( num_exo_blocks / 2)
{
  blocks        = new ElemConnectivity::BlockInfo[max_block_id];
  first_block_id = new int[max_block_id+1];
  
  // The first entry in the block list is unused.
  // The list is indexed by block id, and the first
  // block id is 1.
  blocks[0].elem_count = 0;
  blocks[0].elem_size  = 0;
  blocks[0].elem_type  = INVALID_ELEMENT_TYPE;
  blocks[0].base_type  = INVALID_ELEMENT_TYPE;
  blocks[0].elem_conn  = 0;
  first_block_id[0]    = 0;
  first_block_id[1]    = 1;
  
  int elem_count, elem_size ;
  char elem_type[MAX_STR_LENGTH+1];
  int attr, error;
  ElemConnectivity::BlockInfo* block = blocks + 1;
  for( int i = 1; i < max_block_id; i++ )
  {
    error = ex_get_elem_block( file_id, i, elem_type, 
                               &elem_count, &elem_size, 
                               &attr );
    if( error )
    {
      PRINT_ERROR("Problem reading element block %d.\n",i);
      exit(2);
    }
    
    block->elem_count  = elem_count;
    block->elem_size   = elem_size;
    block->base_type   = get_base_elem_type( elem_type );
    block->elem_type   = get_element_type( block->base_type, elem_type );
    block->elem_conn   = 0;
    first_block_id[i+1]= first_block_id[i] + elem_count;
    block++;
  }
}

ElemConnectivity::~ElemConnectivity()
{
  for( int i = 1; i < max_block_id; i++ )
    delete [] blocks[i].elem_conn;
  delete [] blocks;
  delete [] first_block_id;
}

int ElemConnectivity::get_block_id( int elem_id ) 
{
  int block_id = last_block;
  if( elem_id < first_block_id[last_block] )
    block_id = 0;
  
  while( (block_id < max_block_id) && 
         (first_block_id[block_id+1] <= elem_id) )
    block_id++;
      

  if( block_id > max_block_id )
    block_id = 0;
  if( block_id == 0 )
    return 0;
    
   assert( (elem_id >= first_block_id[block_id]) &&
           (elem_id <  first_block_id[block_id+1]) ); 
  return (last_block = block_id);
}

const int* ElemConnectivity::get_elem_conn( int elem_id, int& block_id )
{
  block_id = get_block_id( elem_id );
  if( block_id < 1 ) return 0;
  ElemConnectivity::BlockInfo& block = blocks[block_id];
  int* elem_conn   = block.elem_conn;
  int  elem_size   = block.elem_size;
  int  elem_count  = block.elem_count;
  int  first_id    = first_block_id[block_id];
  ElemType type    = block.elem_type;
  if( ! elem_conn )
  {
    elem_conn = new int[elem_count*elem_size];
    int error = ex_get_elem_conn( file_id, block_id, elem_conn );
    if( error )
    {
      PRINT_ERROR("Problems reading connectivity info for block %d.\n",block_id);
      exit(2);
    }
    block.elem_conn = elem_conn;
  }
  
  return elem_conn + ( (elem_id - first_id) * elem_size );
}

int ElemConnectivity::get_elem( int node_ids[8], int elem_id, int side_id )
{
  int block_id;
  const int* elem_conn = get_elem_conn( elem_id, block_id );
  if( ! elem_conn ) return -1;
  
  switch( blocks[block_id].base_type )
  {
    case HEXSHELL:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 8 * sizeof(int) );
          return 8;
        case 1:
          memcpy( node_ids, elem_conn + 8, 4 * sizeof(int ) );
          return 4;
        case 2:
          node_ids[0] = elem_conn[11];
          node_ids[1] = elem_conn[10];
          node_ids[2] = elem_conn[9];
          node_ids[3] = elem_conn[8];
          return 4;
        default:
          //fall through to normal hex handler.
          side_id -= 2;
      }
    
    case HEX:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 8 * sizeof(int) );
          return 8;
        case 1:
        case 2:
        case 3:
        case 4:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%4];
          node_ids[2] = elem_conn[side_id%4+4];
          node_ids[3] = elem_conn[side_id+3];
          return 4;
        case 5:
          node_ids[0] = elem_conn[3];
          node_ids[1] = elem_conn[2];
          node_ids[2] = elem_conn[1];
          node_ids[3] = elem_conn[0];
          return 4;
        case 6:
          node_ids[0] = elem_conn[4];
          node_ids[1] = elem_conn[5];
          node_ids[2] = elem_conn[6];
          node_ids[3] = elem_conn[7];
          return 4;
      }
      break;
    
    case TETRA: 
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 4 * sizeof(int) );
          return 4;
        case 1:
        case 2:
        case 3:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%3];
          node_ids[2] = elem_conn[3];
          return 3;
        case 4:
          node_ids[0] = elem_conn[2];
          node_ids[1] = elem_conn[1];
          node_ids[2] = elem_conn[0];
          return 3;
      }
      break;
    
    case PYRAMID:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 5 * sizeof(int) );
          return 5;
        case 1:
        case 2:
        case 3:
        case 4:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%4];
          node_ids[2] = elem_conn[4];
          return 3;
        case 5:
          node_ids[0] = elem_conn[3];
          node_ids[1] = elem_conn[2];
          node_ids[2] = elem_conn[1];
          node_ids[3] = elem_conn[0];
          return 4;
      }
      break;
    
    case WEDGE:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 6 * sizeof(int) );
          return 6;
        case 1:
        case 2:
        case 3:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%3];
          node_ids[2] = elem_conn[side_id%3+3];
          node_ids[3] = elem_conn[side_id+2];
          return 4;
        case 4:
          node_ids[0] = elem_conn[2];
          node_ids[1] = elem_conn[1];
          node_ids[2] = elem_conn[0];
          return 3;
        case 5:
          node_ids[0] = elem_conn[3];
          node_ids[1] = elem_conn[4];
          node_ids[2] = elem_conn[5];
          return 3;
      }
      break;
    
    case TRI:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 3 * sizeof(int) );
          return 3;
        case 1:
        case 2:
        case 3:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%3];
          return 2;
      }
      break;        
    
    case SHELL:
      switch( side_id )
      {
        case 0:
        case 1:
          memcpy( node_ids, elem_conn, 4 * sizeof(int) );
          return 4;
        case 2:
          node_ids[3] = elem_conn[0];
          node_ids[2] = elem_conn[1];
          node_ids[1] = elem_conn[2];
          node_ids[0] = elem_conn[3];
          return 4;
        case 3:
        case 4:
        case 5:
        case 6:
          node_ids[0] = elem_conn[side_id-3];
          node_ids[1] = elem_conn[(side_id-2)%4];
          return 2;
      }
      break;
    
    case QUAD:
      switch( side_id )
      {
        case 0:
          memcpy( node_ids, elem_conn, 4 * sizeof(int) );
          return 4;
        case 1:
        case 2:
        case 3:
        case 4:
          node_ids[0] = elem_conn[side_id-1];
          node_ids[1] = elem_conn[side_id%4];
          return 2;
      }
      break;
    
    case BAR  :
    case BEAM :
    case TRUSS:
      switch( side_id )
      {
        case 0:
          node_ids[0] = elem_conn[0];
          node_ids[1] = elem_conn[1];
          return 2;
        case 1:
          node_ids[0] = elem_conn[0];
          return 1;
        case 2:
          node_ids[0] = elem_conn[1];
          return 1;
      }
      break;
  }
  
  PRINT_ERROR("Invalid element type (\"%s\") or element face id (%d).\n",
    element_type_names[blocks[block_id].elem_type], side_id );
  return 0;  
}

ElemType ElemConnectivity::get_base_elem_type( const char* elem_name )
{
  ElemType result = INVALID_ELEMENT_TYPE;
  if(! strncmp( elem_name, "BAR", 3 ) )
    result = BAR;
  else if(! strncmp( elem_name, "BEAM", 4 ) )
    result = BEAM;
  else if(! strncmp( elem_name, "TRUSS", 5 ) )
    result = TRUSS;
  else if(! strncmp( elem_name, "TRI", 3 ) )
    result = TRI;
  else if(! strncmp( elem_name, "QUAD", 4 ) )
    result = QUAD;
  else if(! strncmp( elem_name, "SHEL", 4 ) )
    result = SHELL;
  else if(! strncmp( elem_name, "TETRA", 5 ) )
    result = TETRA;
  else if(! strncmp( elem_name, "PYRAMID", 7 ) )
    result = PYRAMID;
  else if(! strncmp( elem_name, "WEDGE", 5 ) )
    result = WEDGE;
  else if(! strncmp( elem_name, "HEXSHELL", 8 ) )
    result = HEXSHELL;
  else if(! strncmp( elem_name, "HEX", 3 ) )
    result = HEX;
  
  return result;
}



ElemType ElemConnectivity::get_element_type( const char* elem_name )
{
  return get_element_type( BAR, elem_name );
}
ElemType ElemConnectivity::get_element_type( ElemType base, const char* name )
{
  int result;
  for( result = base; result < INVALID_ELEMENT_TYPE; result++ )
  {
    if( !strcmp( element_type_names[result], name ) )
    {
      break;
    }
  }
  return (ElemType)result;
}

void ElemConnectivity::free_memory()
{
  for( int i = 1; i < max_block_id; i++ )
  {
    //if( i == last_block ) continue;
    ElemConnectivity::BlockInfo& block = blocks[i];
    if( block.elem_conn )
    {
      delete [] block.elem_conn;
      block.elem_conn = 0;
    }
  }
}

int ElemConnectivity::elem_dimension( ElemType type )
{
  if( type < BAR )
  {
    return 0;
  }
  else if( type < TRI )
  {
    return 1;
  }
  else if( type < TETRA )
  {
    return 2;
  }
  else if( type < INVALID_ELEMENT_TYPE )
  {
    return 3;
  }
  else
  {
    return 0;
  }
}
