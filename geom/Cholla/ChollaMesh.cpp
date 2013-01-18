
#include "ChollaMesh.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"

ChollaMesh::ChollaMesh()
{
}

ChollaMesh::~ChollaMesh()
{
}

FacetEntity* ChollaMesh::convert_to_facet(SMD::SMDElement elem)
{
  return reinterpret_cast<FacetEntity*>(elem.mHandle);
}

CubitPoint* ChollaMesh::convert_to_point(SMD::SMDNode node)
{
  return reinterpret_cast<CubitPoint*>(node.mHandle);
}

SMD::SMDElement ChollaMesh::convert_from_facet(FacetEntity* f)
{
  SMD::SMDElement e;
  e.mHandle = reinterpret_cast<size_t>(f);
  return e;
}

SMD::SMDNode ChollaMesh::convert_from_point(CubitPoint* p)
{
  SMD::SMDNode n;
  n.mHandle = reinterpret_cast<size_t>(p);
  return n;
}


SMD::ErrorCode ChollaMesh::get_ids(size_t num, const SMD::SMDNode nodes[], int ids[])
{
  for(size_t i=0; i<num; i++)
  {
    ids[i] = convert_to_point(nodes[i])->id();
  }
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_ids(size_t num, const SMD::SMDElement elems[], int ids[])
{
  CubitFacet* tri;
  CubitFacetEdge* edge;

  for(size_t i=0; i<num; i++)
  {
    FacetEntity* f = convert_to_facet(elems[i]);
    if((tri = dynamic_cast<CubitFacet*>(f)))
    {
      ids[i] = tri->id();
    }
    else if((edge = dynamic_cast<CubitFacetEdge*>(f)))
    {
      ids[i] = edge->id();
    }
    else
    {
      ids[i] = 0;
    }
  }
  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::get_node_handle(int id, SMD::SMDNode &node)
{
  return SMD::ERR_NOT_IMPLEMENTED;
}

SMD::ErrorCode ChollaMesh::get_element_handle(SMD::ElementType type, int id, SMD::SMDElement &node)
{
  return SMD::ERR_NOT_IMPLEMENTED;
}

SMD::ErrorCode ChollaMesh::get_element_types(int num_elems, const SMD::SMDElement elem_handle_array[], 
                                 SMD::ElementType type_array[], bool*)
{
  for(int i=0; i<num_elems; i++)
  {
    FacetEntity* f = convert_to_facet(elem_handle_array[i]);
    if(dynamic_cast<CubitFacet*>(f))
    {
      type_array[i] = SMD::TRI;
    }
    else if(dynamic_cast<CubitFacetEdge*>(f))
    {
      type_array[i] = SMD::EDGE;
    }
    else
    {
      type_array[i] = SMD::NO_ELEMENT_TYPE;
    }
  }
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_element_dimensions(int num_elems, const SMD::SMDElement elem_handle_array[],
                                      int dimensions[])
{
  for(int i=0; i<num_elems; i++)
  {
    FacetEntity* f = convert_to_facet(elem_handle_array[i]);
    if(dynamic_cast<CubitFacet*>(f))
    {
      dimensions[i] = 2;
    }
    else if(dynamic_cast<CubitFacetEdge*>(f))
    {
      dimensions[i] = 1;
    }
    else
    {
      dimensions[i] = 0;
    }
  }
  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::get_element_length(SMD::SMDElement elem_handle, double& length)
{
  CubitFacetEdge* f = dynamic_cast<CubitFacetEdge*>(convert_to_facet(elem_handle));
  if(!f)
  {
    return SMD::ERR_INVALID_ENTITY_HANDLE;
  }

  length = f->length();
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_element_area(SMD::SMDElement elem_handle, double& area)
{
  CubitFacet* f = dynamic_cast<CubitFacet*>(convert_to_facet(elem_handle));
  if(!f)
  {
    return SMD::ERR_INVALID_ENTITY_HANDLE;
  }

  area = f->area();
  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::get_element_normal(SMD::SMDElement elem_handle, double normal_vector[3])
{
  CubitFacet* f = dynamic_cast<CubitFacet*>(convert_to_facet(elem_handle));
  if(!f)
  {
    return SMD::ERR_INVALID_ENTITY_HANDLE;
  }

  CubitVector n = f->normal();
  normal_vector[0] = n.x();
  normal_vector[1] = n.y();
  normal_vector[2] = n.z();
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_element_centroid(SMD::SMDElement elem_handle, double centroid[3])
{
  CubitFacet* f = dynamic_cast<CubitFacet*>(convert_to_facet(elem_handle));
  if(!f)
  {
    return SMD::ERR_INVALID_ENTITY_HANDLE;
  }

  CubitVector c = f->center();
  centroid[0] = c.x();
  centroid[1] = c.y();
  centroid[2] = c.z();
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], float coords[][3])
{
  for(int i=0; i<num_nodes; i++)
  {
    CubitPoint* p = convert_to_point(node_handle_array[i]);
    coords[i][0] = p->x();
    coords[i][1] = p->y();
    coords[i][2] = p->z();
  }
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], double coords[][3])
{
  for(int i=0; i<num_nodes; i++)
  {
    CubitPoint* p = convert_to_point(node_handle_array[i]);
    coords[i][0] = p->x();
    coords[i][1] = p->y();
    coords[i][2] = p->z();
  }
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_coords(int num_nodes, const SMD::SMDNode node_handle_array[], 
                          SMD::Real x_coords[], SMD::Real y_coords[], SMD::Real z_coords[])
{
  for(int i=0; i<num_nodes; i++)
  {
    CubitPoint* p = convert_to_point(node_handle_array[i]);
    x_coords[i] = p->x();
    y_coords[i] = p->y();
    z_coords[i] = p->z();
  }
  return SMD::STATUS_OK;
}

SMD::ErrorCode ChollaMesh::get_connectivity(
  int num_elems,
  const SMD::SMDElement elements[],
  int& node_array_size,
  SMD::SMDNode node_handles[])
{
  CubitFacet* tri;
  CubitFacetEdge* edge;

  int offset = 0;
  for(int i=0; i<num_elems; i++)
  {
    FacetEntity* f = convert_to_facet(elements[i]);
    if((tri = dynamic_cast<CubitFacet*>(f)))
    {
      if(offset + 3 > node_array_size)
      {
        return SMD::ERR_ARRAY_TOO_SMALL;
      }

      node_handles[offset++] = convert_from_point(tri->point(0));
      node_handles[offset++] = convert_from_point(tri->point(1));
      node_handles[offset++] = convert_from_point(tri->point(2));
    }
    else if((edge = dynamic_cast<CubitFacetEdge*>(f)))
    {
      if(offset + 2 > node_array_size)
      {
        return SMD::ERR_ARRAY_TOO_SMALL;
      }
      node_handles[offset++] = convert_from_point(edge->point(0));
      node_handles[offset++] = convert_from_point(edge->point(1));
    }
    else
    {
      return SMD::ERR_INVALID_ENTITY_HANDLE;
    }
  }
  node_array_size = offset;
  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::get_expanded_connectivity(
  int num_elems,
  const SMD::SMDElement elements[],
  unsigned int nodes_per_elem,
  int node_array_size,
  SMD::SMDNode node_handles[])
{
  return SMD::ERR_NOT_IMPLEMENTED;
}

SMD::ErrorCode ChollaMesh::get_expanded_connectivity(
  SMD::SMDElement element,
  unsigned int& nodes_per_elem,
  int node_array_size,
  SMD::SMDNode node_handles[])
{
  SMD::ErrorCode ret = get_connectivity(1, &element, node_array_size, node_handles);
  nodes_per_elem = node_array_size;
  return ret;
}

SMD::ErrorCode ChollaMesh::create_nodes(
  unsigned int num_nodes,
  const SMD::Real coords[][3],
  SMD::SMeshOwner owner,
  SMD::SMDNode node_handles[])
{
  for(unsigned int i=0; i<num_nodes; i++)
  {
    node_handles[i] = convert_from_point(new CubitPointData(coords[i][0], coords[i][1], coords[i][2]));
  }
  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::create_elements(
  SMD::ElementType type,
  unsigned int num_nodes,
  const SMD::SMDNode node_handle_array[],
  unsigned int num_elements,
  unsigned int num_nodes_per_element,
  const unsigned int connectivity[],
  SMD::SMeshOwner new_entities_owner,
  SMD::SMDElement *created_element_handles)
{
  if(type == SMD::EDGE)
  {
    unsigned int conn_offset=0;
    for(unsigned int i=0; i<num_elements; i++, conn_offset += num_nodes_per_element)
    {
      CubitPoint* pts[2];
      pts[0] = convert_to_point(node_handle_array[connectivity[conn_offset]]);
      pts[1] = convert_to_point(node_handle_array[connectivity[conn_offset+1]]);
      created_element_handles[i] = convert_from_facet(new CubitFacetEdgeData(pts[0], pts[1]));
    }
    return SMD::STATUS_OK;
  }
  else if(type == SMD::TRI)
  {
    unsigned int conn_offset=0;
    for(unsigned int i=0; i<num_elements; i++, conn_offset += num_nodes_per_element)
    {
      CubitPoint* pts[2];
      pts[0] = convert_to_point(node_handle_array[connectivity[conn_offset]]);
      pts[1] = convert_to_point(node_handle_array[connectivity[conn_offset+1]]);
      pts[2] = convert_to_point(node_handle_array[connectivity[conn_offset+2]]);
      created_element_handles[i] = convert_from_facet(new CubitFacetData(pts[0], pts[1], pts[2]));
    }
    return SMD::STATUS_OK;
  }

  return SMD::ERR_INVALID_ENTITY_TYPE;
}


SMD::ErrorCode ChollaMesh::delete_node(SMD::SMDNode node)
{
  CubitPoint* f = convert_to_point(node);
  if(!f)
    return SMD::ERR_INVALID_ENTITY_HANDLE;

  delete f;

  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::delete_element(SMD::SMDElement element)
{
  FacetEntity* f = convert_to_facet(element);
  if(!f)
    return SMD::ERR_INVALID_ENTITY_HANDLE;

  delete f;

  return SMD::STATUS_OK;
}


SMD::ErrorCode ChollaMesh::reverse_element_connectivity(SMD::SMDElement element)
{
  FacetEntity* f = convert_to_facet(element);
  CubitFacet* tri;
  CubitFacetEdge* edge;

  if((tri = dynamic_cast<CubitFacet*>(f)))
  {
    tri->flip();
    return SMD::STATUS_OK;
  }
  else if((edge = dynamic_cast<CubitFacetEdge*>(f)))
  {
    edge->flip();
    return SMD::STATUS_OK;
  }
  return SMD::ERR_INVALID_ENTITY_HANDLE;
}

SMD::ErrorCode ChollaMesh::find_element(SMD::ElementType type, unsigned int num_pts, const SMD::SMDNode nodes[], 
                                SMD::SMDElement& found_element)
{
  if(type == SMD::EDGE && num_pts == 2)
  {
    CubitFacetEdge* edge = convert_to_point(nodes[0])->shared_edge(convert_to_point(nodes[1]));
    found_element = convert_from_facet(edge);
    return SMD::STATUS_OK;
  }
  else if(type == SMD::TRI && num_pts == 3)
  {
    DLIList<CubitFacet*> facets[3];
    convert_to_point(nodes[0])->facets(facets[0]);
    convert_to_point(nodes[1])->facets(facets[1]);
    convert_to_point(nodes[2])->facets(facets[2]);
    facets[0].intersect_unordered(facets[1]);
    facets[0].intersect_unordered(facets[2]);

    if(facets[0].size())
    {
      found_element = convert_from_facet(facets[0][0]);
    }
    return SMD::STATUS_OK;
  }
  return SMD::ERR_INVALID_ENTITY_TYPE;
}

