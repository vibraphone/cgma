
#ifndef ChollaMesh_hpp
#define ChollaMesh_hpp

#include "SMDInterface.hpp"

class FacetEntity;
class CubitPoint;

// a SMD interface to Cholla data structures
class ChollaMesh : public SMD::MeshData
{
  public:
    ChollaMesh();
    ~ChollaMesh();

    FacetEntity* convert_to_facet(SMD::SMDElement elem);
    CubitPoint* convert_to_point(SMD::SMDNode node);
    SMD::SMDElement convert_from_facet(FacetEntity* f);
    SMD::SMDNode convert_from_point(CubitPoint* p);


    SMD::ErrorCode get_ids(size_t num, const SMD::SMDNode nodes[], int ids[]);
    SMD::ErrorCode get_ids(size_t num, const SMD::SMDElement elems[], int ids[]);
    
    SMD::ErrorCode get_node_handle(int id, SMD::SMDNode &node);
    SMD::ErrorCode get_element_handle(SMD::ElementType type, int id, SMD::SMDElement &elem);
    SMD::ErrorCode get_element_types(int num_elems, const SMD::SMDElement elem_handle_array[], 
                                     SMD::ElementType type_array[], bool* same_type);
    SMD::ErrorCode get_element_dimensions(int num_elems, const SMD::SMDElement elem_handle_array[],
                                          int dimensions[]);

    SMD::ErrorCode get_element_length(SMD::SMDElement elem_handle, double& length);
    SMD::ErrorCode get_element_area(SMD::SMDElement elem_handle, double& area);
    
    SMD::ErrorCode get_element_normal(SMD::SMDElement elem_handle, double normal_vector[3]);
    
    SMD::ErrorCode get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], float coords[][3]);
    SMD::ErrorCode get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], double coords[][3]);
    SMD::ErrorCode get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], 
                              SMD::Real x_coords[], SMD::Real y_coords[], SMD::Real z_coords[]);
    SMD::ErrorCode get_connectivity(
      int num_elems,
      const SMD::SMDElement elements[],
      int& node_array_size,
      SMD::SMDNode node_handles[]);
    
    SMD::ErrorCode get_expanded_connectivity(
      int num_elems,
      const SMD::SMDElement elements[],
      unsigned int nodes_per_elem,
      int node_array_size,
      SMD::SMDNode node_handles[]);
    
    SMD::ErrorCode get_expanded_connectivity(
      SMD::SMDElement elements,
      unsigned int& nodes_per_elem,
      int node_array_size,
      SMD::SMDNode node_handles[]);
    
    SMD::ErrorCode create_nodes(
      unsigned int num_nodes,
      const SMD::Real coords[][3],
      SMD::SMeshOwner owner,
      SMD::SMDNode node_handles[]);
    
    SMD::ErrorCode create_elements(
      SMD::ElementType type,
      unsigned int num_nodes,
      const SMD::SMDNode node_handle_array[],
      unsigned int num_elements,
      unsigned int num_nodes_per_element,
      const unsigned int connectivity[],
      SMD::SMeshOwner new_entities_owner,
      SMD::SMDElement *created_element_handles);
    
    SMD::ErrorCode delete_node(SMD::SMDNode node);

    SMD::ErrorCode delete_element(SMD::SMDElement element);
    
    SMD::ErrorCode reverse_element_connectivity(SMD::SMDElement element);
  
    SMD::ErrorCode find_element(SMD::ElementType type, unsigned int num_pts, const SMD::SMDNode nodes[], 
                                SMD::SMDElement& found_element);

};

#endif // ChollaMesh_hpp

