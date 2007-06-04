#ifndef ELEM_CONNECTIVITY_HPP
#define ELEM_CONNECTIVITY_HPP

enum ElemType {
  BAR=0, BAR2,   BAR3,
  BEAM,  BEAM2,  BEAM3,
  TRUSS, TRUSS2, TRUSS3,
  TRI,   TRI3,   TRI4,   TRI6,   TRI7,
  SHELL, SHELL4, SHELL5, SHELL8, SHELL9,
  QUAD,  QUAD4,  QUAD5,  QUAD8,  QUAD9,
  TETRA,  TETRA4,  TETRA5,  TETRA8,   TETRA10,  TETRA14,  TETRA15,
  PYRAMID,PYRAMID5,PYRAMID6,PYRAMID10,PYRAMID13,PYRAMID18,PYRAMID19,
  WEDGE,  WEDGE6,  WEDGE7,  WEDGE11,  WEDGE15,  WEDGE20,  WEDGE21,
  HEX,    HEX8,    HEX9,    HEX14,    HEX20,    HEX26,    HEX27,
  HEXSHELL, INVALID_ELEMENT_TYPE 
};


class ElemConnectivity
{
    struct BlockInfo
    {
      int  elem_count;  // number of elements in the block
      int  elem_size;   // number of nodes per element
      ElemType elem_type;  // type of element in block
      ElemType base_type;  // base elem type
      int* elem_conn;   // element connectivity array
    };
  
  public:
    
    ElemConnectivity( int exo_file_id, int num_exo_blocks );
    
    ~ElemConnectivity();
    
    int get_elem( int node_ids[8], int elem_id, int side_id = 0 );
    //- Passes back the node_ids for the specified element.  
    //- elem_id is the id of the desired element.
    //- side_id, if non-zero, specifies a side of the element.
    //-  For example, specifying a hex elem_id and a non-zero
    //-  side id will get the specified face of the hex.  A non-
    //-  zero side_id for a quad will get an edge, etc.
    //- The return value is the number of nodes passed back in
    //- node_ids.  If an error is encountered, zero is returned.
    //-
    //- This function disregards all mid nodes.  For example, for
    //- a HEX27 element, only the 8 corner nodes are returned.
    
    int get_block_id( int elem_id ) ;
    //- Get the id of the block containing th specified element.
    
    static ElemType get_base_elem_type( const char* elem_name );
    //- Get the base type (e.g. HEX for HEX27).
    
    static ElemType get_element_type( const char* elem_name );
    static ElemType get_element_type( ElemType base_type,
                                      const char* elem_name );
    //- Get enum for passed string element name.
    
    static const char* const element_type_names[];
    //- Array of element names, indexed by ElemType.
    
    int first_elem_id( int block_id ) const  
      { return first_block_id[block_id]; }
    //- Get the Id of the first element in a block.
    
    int elem_count( int block_id ) const 
      { return blocks[block_id].elem_count; }
    //- Get the number of elements in a block.
    
    ElemType elem_type( int block_id ) const 
      { return blocks[block_id].elem_type; }
    //- Get the type of element in a block
    
    ElemType base_type( int block_id ) const
      { return blocks[block_id].base_type; }
    
    const int block_count;  
      // number of element blocks in the file
    
    const int max_block_id; 
      // block_count + 1
    
    const int file_id;
      // Exodus file id.
      
    void free_memory();
    
    static int elem_dimension( ElemType type );

  private:
    
    const int* get_elem_conn( int elem_id, int& block_id );
      // Get a pointer to the connectivity information for
      // the specified element.

    BlockInfo* blocks;     
      // Metadata for each block.
      // Array is indexed by block id, so first entry is empty.
      
    int* first_block_id;
      // The id fo the first element in each block.
      // The array is indexed by block id.  The first
      // and last entries do not correspond to any
      // block, and contain the values 0 and max_elem_id + 1, 
      // respectively.
      
    int last_block;
      // block the last requested element was in.
};

#endif

