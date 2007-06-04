#include "exodusII.h"
#include "CGMApp.hpp"
#include "AppUtil.hpp"
#include "GeometryQueryTool.hpp"
#include "FacetModifyEngine.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "FacetSurface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "ElemConnectivity.hpp"
#include <vector>

#ifndef CANT_USE_STD
using std::vector;
#endif

const char usage[] = 
  "Usage: %s <exodus_file> -body {-block <id_list> | -sideset <id_list>}\n"
  "                   [-body {-block <id_list> | -sideset <id_list}]... \n"
  "                   [-angle <degrees>]\n";

extern "C" { void gl_cleanup() { } }

  
class ExoLoader
{
  public:
  
  ExoLoader( int argc, char* argv[] );
  
  ~ExoLoader();
  
  private:
  
  void do_interface();
  
  bool open_exodus_file( const char* file_name );
  
  CubitPointData* get_point( int point_id )
  {
    --point_id;
    CubitPointData** ptr = point_list + point_id;
    if( !*ptr ) 
      *ptr = new CubitPointData( x_coords[point_id],
                                 y_coords[point_id],
                                 z_coords[point_id] );
    return *ptr;
  }
  
  CubitFacet* find_tri ( CubitPoint*, CubitPoint*, CubitPoint*);
  void add_tri ( CubitPoint*, CubitPoint*, CubitPoint*);
  void add_quad( CubitPoint*, CubitPoint*, CubitPoint*, CubitPoint*);
  
  bool load_block(   int block_id   );
  bool load_sideset( int sideset_id );
  
  bool make_curves(   vector<int>& edge_set );
  bool make_surfaces( vector<int>& face_set );
  bool create_geometry();
  
  static bool incomplete_compare( const char* whole, 
                                  const char* part, 
                                  int min_match = 2);
                                  
  static void eval_surface( Surface* surf_ptr );
  

  
  // Node coordinates
  double *x_coords;
  double *y_coords;
  double *z_coords;
  CubitPointData** point_list;
  
  // ExodusII file Id
  int exo_file_id;
  
  // IDs of active entities
  vector<int> active_block_list;
  vector<int> active_sides_list;
  
  // Exodus elements
  ElemConnectivity* elements;
  
  // Exodus info.
  int exo_node_count;
  int exo_elem_count;
  int exo_elem_block_count;
  int exo_node_set_count;
  int exo_side_set_count;
  
  // CGM interface
  GeometryQueryTool* geom_tool;
  FacetModifyEngine* facet_tool;
  double angle;
  
  // Geometry
  vector<RefEntity*> geometry;
  DLIList<CubitFacet*> tri_facets;
//  DLIList<CubitQuadFacet*> quad_facets;
  
  // print debug output
  int debug;
};



int main( int argc, char* argv[] )
{
  delete new ExoLoader( argc, argv );
  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during exoIIfacet session.\n");
  }
  return ret_val;
}

ExoLoader::ExoLoader( int argc, char* argv[] )
{
  // Initialize some stuff
  x_coords = y_coords = z_coords = 0;
  exo_file_id = -1;
  elements = 0;
  angle = 140.;
  debug = 0;
  
  // Intialize CGM
  AppUtil::instance()->startup( argc, argv );
  CGMApp::instance()->startup( argc, argv );
  facet_tool = FacetModifyEngine::instance();
  
  // Parse arguments
  
  if( argc < 2 )
  {
    PRINT_ERROR(usage,argv[0]);
    return;
  }
  
  if( !open_exodus_file( argv[1] ) )
  {
    PRINT_ERROR(usage,argv[0]);
    return;
  }
  
  
  int i = 2;
  while( i < argc )
  {
    if( incomplete_compare( "-block", argv[i] ) )
    {
      i++;
      while( i < argc )
      {
        char* ptr = argv[i];
        long block_id = strtol(ptr, &ptr, 0);
        if( (ptr != argv[i]) && (*ptr == '\0') && 
            (block_id > 0) && (block_id < INT_MAX) )
        {
          i++;
          load_block( (int)block_id );
        }
        else
        {
          break;
        }
      }
    }
    else if( incomplete_compare( "-sideset", argv[i] ) )
    {
      i++;
      while( i < argc )
      {
        char* ptr = argv[i];
        long block_id = strtol(ptr, &ptr, 0);
        if( (ptr != argv[i]) && (*ptr == '\0') && 
            (block_id > 0) && (block_id < INT_MAX) )
        {
          i++;
          load_sideset( (int)block_id );
        }
        else
        {
          break;
        }
      }
    }
    else if( incomplete_compare( "-angle", argv[i] ) )
    {
      if( (i+1) < argc )
      {
        i++;
        char* ptr = argv[i];
        angle = strtod( ptr, &ptr );
        if( (ptr != argv[i]) && (*ptr == '\0') )
        {
          i++;
        }
        else
        {
          angle = 140;
        }
      }
    }
    else if( incomplete_compare( "-body", argv[i] ) )
    {
      i++;
      create_geometry();
    }
    else if( incomplete_compare( "-debug", argv[i] ) )
    {
      debug = 1;
      i++;
    }
    else
    {
      PRINT_ERROR("Unrecognized option \"%s\".\n",argv[i] );
      PRINT_ERROR(usage,argv[0]);
      return;
    }
  }
  
  // Do whatever the interface is going to be.
  create_geometry();
}

ExoLoader::~ExoLoader()
{
  delete [] x_coords;
  delete [] y_coords;
  delete [] z_coords;
  delete [] point_list;
  delete elements;
}
  

bool ExoLoader::open_exodus_file( const char* filename )
{
  // open the file
  assert( exo_file_id == -1 );
  assert( !elements );
  
  int exo_cpu_real_size = sizeof(double);
  int exo_io_real_size = 0;
  float exo_version;
  exo_file_id = ex_open( filename, EX_READ,
                         &exo_cpu_real_size,
                         &exo_io_real_size,
                         &exo_version );
  if( exo_file_id < 0 )
  {
    PRINT_ERROR("Could not read file \"%s\"\n",filename);
    return false;
  }
  
  char exo_file_title[MAX_LINE_LENGTH];
  int exo_dimension;
  int exo_error = ex_get_init( exo_file_id, 
                               exo_file_title, 
                               &exo_dimension,
                               &exo_node_count, 
                               &exo_elem_count,
                               &exo_elem_block_count,
                               &exo_node_set_count,
                               &exo_side_set_count );
  
  if( exo_error )
  {
    PRINT_ERROR("Problems reading file \"%s\"\n",filename);
    return false;
  }
  
  if( debug )
  {
    PRINT_INFO("Exodus File:\n"
               "\t id     = %d\n"
               "\t title  = %s\n"
               "\t dim    = %d\n"
               "\t nodes  = %d\n"
               "\t elems  = %d\n"
               "\t blocks = %d\n"
               "\t sides  = %d\n",
               exo_file_id,
               exo_file_title,
               exo_dimension,
               exo_node_count,
               exo_elem_count,
               exo_elem_block_count,
               exo_side_set_count);
  }
  
  if( exo_dimension != 3 )
  {
    PRINT_WARNING("Exodus data is only %d dimensional.  Assuming zero for"
                  " coordinate values for other dimensions.\n", exo_dimension );
  }
  
  x_coords = new double[exo_node_count];
  y_coords = new double[exo_node_count];
  z_coords = new double[exo_node_count];
  if( exo_dimension < 2 ) memset( y_coords, 0, exo_node_count * sizeof(double) );
  if( exo_dimension < 3 ) memset( z_coords, 0, exo_node_count * sizeof(double) );
  exo_error = ex_get_coord( exo_file_id, x_coords, y_coords, z_coords );
  if( exo_error )
  {
    PRINT_ERROR("Problems reading node coordinates.\n");
    delete [] x_coords;
    delete [] y_coords;
    delete [] z_coords;
    x_coords = y_coords = z_coords = 0;
    ex_close( exo_file_id );
    exo_file_id = -1;
    exo_node_count = exo_elem_count = exo_elem_block_count = 0;
    return false;
  }
  else
  {
    PRINT_INFO("Read %d nodes from \"%s\".\n", exo_node_count, filename );
  }
  
  point_list = new CubitPointData*[exo_node_count];
  memset( point_list, 0, exo_node_count * sizeof( void* ) );
  
  elements = new ElemConnectivity( exo_file_id, exo_elem_block_count );
  return true;
}

bool ExoLoader::incomplete_compare( const char* whole, const char* part, int min )
{
  const char* ptr1 = whole;
  const char* ptr2 = part;
  while( *ptr1 && (*ptr1 == *ptr2) )
  {
    ptr1++;
    ptr2++;
  }
  
  return *ptr2 || (ptr2 - part) < min ? false : true;
}

  
bool ExoLoader::load_block( int block_id )
{
  if( (block_id < 1) || (block_id > exo_elem_block_count) )
  {
    PRINT_ERROR("Invalid block Id %d\n",block_id);
    return false;
  }
  
  ElemType type = elements->base_type( block_id );
  
  int elem_count = elements->elem_count(block_id);
  int id_offset  = elements->first_elem_id(block_id);
  vector<int> elem_list(elem_count * 2);
  int face_id = type == HEXSHELL ? 1 : 0;
  for( int i = 0; i < elem_count; i++ )
  {
    elem_list[i*2] = id_offset++;
    elem_list[i*2+1] = face_id;
  }
  
  if( debug )
    PRINT_INFO("Block %d: \n\t type  = %s\n\t elems = %d\n\t first = %d\n",
      block_id, ElemConnectivity::element_type_names[type], id_offset );
  
  switch( type )
  {
    case BAR: case BEAM: case TRUSS:
      return make_curves( elem_list );
    case TRI: case SHELL: case QUAD: case HEXSHELL:
      return make_surfaces( elem_list );
    default:
      PRINT_ERROR("Don't know how to create geometry from %ss.\n",
        elements->element_type_names[elements->elem_type(block_id)] );
      return false;
  }
}
  
bool ExoLoader::load_sideset( int sideset_id )
{
  if( (sideset_id < 1) || (sideset_id > exo_side_set_count) )
  {
    PRINT_ERROR("Invalid SideSet Id %d\n",sideset_id);
    return false;
  }
  
  int elem_count, junk;
  int exo_error = ex_get_side_set_param( exo_file_id, sideset_id, 
                                         &elem_count, &junk );
  if( exo_error )
  {
    PRINT_ERROR("Could not read sideset %d.\n",sideset_id);
    return false;
  }
  
  int* elem_list = new int[elem_count];
  int* side_list = new int[elem_count];
  exo_error = ex_get_side_set( exo_file_id, sideset_id, elem_list, side_list );
  if( exo_error )
  {
    PRINT_ERROR("Could not read sideset %d.\n",sideset_id);
    delete [] elem_list;
    delete [] side_list;
    return false;
  }
  
  vector<int> one_dim_list(elem_count*2), two_dim_list(elem_count*2);
  int one_d_count = 0;
  int two_d_count = 0;
  for( int i = 0; i < elem_count; i++ )
  {
    int dim = elements->elem_dimension(
                elements->elem_type( 
                  elements->get_block_id( 
                    elem_list[i] ) ) );
    if( side_list[i] != 0 )
      dim--;
  
    if( dim == 1 )
    {
      one_dim_list[one_d_count++] = elem_list[i];
      one_dim_list[one_d_count++] = side_list[i];
    }
    else if( dim == 2 )
    {
      two_dim_list[two_d_count++] = elem_list[i];
      two_dim_list[two_d_count++] = side_list[i];
    }
    else
    {
      PRINT_ERROR("Don't know how to create geometry from %ss.\n"
                  "Aborting sideset geometry creation.\n",
        elements->element_type_names[elements->elem_type(sideset_id)] );
      delete [] elem_list;
      delete [] side_list;
      return false;
    }
  }
  
  one_dim_list.resize( one_d_count );
  two_dim_list.resize( two_d_count );
  
  if( debug )
    PRINT_INFO("Sideset %d: \n"
               "\t total elems = %d\n"
               "\t 1D elems    = %d\n"
               "\t 2D elems    = %d\n",
               sideset_id, elem_count, 
               one_dim_list.size(), 
               two_dim_list.size() );
  
  delete [] elem_list;
  delete [] side_list;

  if( one_dim_list.size() && ! make_curves( one_dim_list ) )
  {
    PRINT_ERROR("Curve creation failed.\n");
    return false;
  }
  
  if( two_dim_list.size() && ! make_surfaces( two_dim_list ) )
  {
    PRINT_ERROR("Surface creation failed.\n");
    return false;
  }
  
  return true;
}


bool ExoLoader::make_surfaces( vector<int>& face_set )
{
  // Create facets
  
  DLIList<Surface*> surface_list;
  
  assert( (face_set.size() % 2) == 0 );
  
  int node_list[8];
  for( int i = 0; i < face_set.size(); i+=2 )
  {
    int elem_id = face_set[i];
    int side_id = face_set[i+1];
    int num_nodes = elements->get_elem( node_list, elem_id, side_id );
    
    if( num_nodes == 3 )
    {
       CubitPointData* pt1 = get_point( node_list[0] );
       CubitPointData* pt2 = get_point( node_list[1] );
       CubitPointData* pt3 = get_point( node_list[2] );
       assert( pt1 && pt2 && pt3 );
       
       add_tri( pt1, pt2, pt3 );
    }
    else if( num_nodes == 4 )
    {
       CubitPointData* pt1 = get_point( node_list[0] );
       CubitPointData* pt2 = get_point( node_list[1] );
       CubitPointData* pt3 = get_point( node_list[2] );
       CubitPointData* pt4 = get_point( node_list[3] );
       assert( pt1 && pt2 && pt3 && pt4 );
       
       //CubitQuadFacetData* facet = new CubitQuadFacetData( pt1, pt2, pt3, pt4 );
       //quad_facets.append( facet );
    
       add_quad( pt1, pt2, pt3, pt4 );
    }
    else if( num_nodes == 0 )
    {
      continue;
    }
    else 
    {
      assert(0);
    }
  }
  
  return true;
}

bool ExoLoader::create_geometry()
{ 
  DLIList<Surface*> surface_list;
  
  if( tri_facets.size() )
  {
    DLIList<CubitPoint*> pt_list;
    DLIList<Surface*> temp_list;
    CubitStatus s = facet_tool->build_facet_surface( NULL, tri_facets, 
      pt_list, angle, 0, CUBIT_FALSE, CUBIT_FALSE, temp_list );
    if( !s )
    {
      PRINT_ERROR( "Surface creation from triangles failed.\n");
    }
    surface_list += temp_list;
    tri_facets.clean_out();
  }
//  if( quad_facets.size() )
//  {
//    DLIList<CubitPoint*> pt_list;
//    DLIList<Surface*> temp_list;
//    CubitStatus s = facet_tool->make_facet_surface( quad_facets,
//      pt_list, angle, 0, CUBIT_FALSE, CUBIT_FALSE, temp_list );
//    if( !s )
//    {
//      PRINT_ERROR( "Surface creation from triangles failed.\n");
//    }
//    surface_list += temp_list;
//    quad_facets.clean_out();
//  }
  
  if( surface_list.size() )
  {
    ShellSM* new_shell = 0;
    CubitStatus s = facet_tool->make_facet_shell( surface_list, new_shell );
    assert( s && new_shell );
    
    DLIList<ShellSM*> shell_list;
    shell_list.append( new_shell );
    Lump* new_lump = 0;
    s = facet_tool->make_facet_lump( shell_list, new_lump );
    assert( s && new_lump );
    
    DLIList<Lump*> lump_list;
    lump_list.append( new_lump );
    BodySM* new_body = 0;
    s = facet_tool->make_facet_body( lump_list, new_body );
    assert( s && new_body );
    
    Body* body = geom_tool->make_Body( new_body );
    assert( body != NULL );
    geometry.push_back( body );
    
    PRINT_INFO("Created body %d containing:\n",body->id());
    
    DLIList<RefFace*> face_list;
    DLIList<RefEdge*> edge_list;
    DLIList<RefVertex*> vtx_list;
    body->ref_faces( face_list );
    body->ref_edges( edge_list );
    body->ref_vertices( vtx_list );
    
    if( face_list.size() )
    {
      PRINT_INFO("\tSurfaces ");
      for( int f = face_list.size(); f > 1 ; f--)
        PRINT_INFO("%d, ", face_list.get_and_step()->id() );
      PRINT_INFO("%d\n",face_list.get_and_step()->id() );
    }
    
    if( edge_list.size() )
    {
      PRINT_INFO("\tCurves ");
      for( int e = edge_list.size(); e > 1 ; e--)
        PRINT_INFO("%d, ", edge_list.get_and_step()->id() );
      PRINT_INFO("%d\n",edge_list.get_and_step()->id() );
    }
    
    if( vtx_list.size() )
    {
      PRINT_INFO("\tVertices ");
      for( int v = vtx_list.size(); v > 1 ; v--)
        PRINT_INFO("%d, ", vtx_list.get_and_step()->id() );
      PRINT_INFO("%d\n",vtx_list.get_and_step()->id() );
    }
    
  }
  
  for( int s = surface_list.size(); s > 0; s-- )
  {
    eval_surface( surface_list.get_and_step() );
  }
  
  return surface_list.size() ? true : false;
}

bool ExoLoader::make_curves( vector<int>& edge_list )
{
  return false;
}

CubitFacet* 
ExoLoader::find_tri( CubitPoint* pt1, CubitPoint* pt2, CubitPoint* pt3 )
{
  DLIList<CubitFacet*> pt_facets;
  pt1->facets( pt_facets );
  for( int i = pt_facets.size(); i--; )
  {
    CubitFacet* facet = pt_facets.get_and_step();
    if( facet->contains( pt2 ) &&
        facet->contains( pt3 ) )
      return facet;
  }
  
  return 0;
}

void 
ExoLoader::add_tri( CubitPoint* pt1, CubitPoint* pt2, CubitPoint* pt3 )
{
  CubitFacet* result = find_tri( pt1, pt2, pt3 );
 
  if( !result )
  {
    result = new CubitFacetData( pt1, pt2, pt3 );
    tri_facets.append( result );
  }
  else 
  {
    tri_facets.append_unique( result );
  }
}

void 
ExoLoader::add_quad( CubitPoint* pt1, CubitPoint* pt2, 
                     CubitPoint* pt3, CubitPoint* pt4 )
{
  CubitFacet* facet1 = find_tri( pt1, pt2, pt3 );
  CubitFacet* facet2 = find_tri( pt3, pt4, pt1 );
  CubitFacet* facet3 = find_tri( pt1, pt2, pt4 );
  CubitFacet* facet4 = find_tri( pt2, pt3, pt4 );
  
  if( facet1 && facet2 )
  {
    tri_facets.append_unique( facet1 );
    tri_facets.append_unique( facet2 );
  }
  else if( facet3 && facet4 )
  {
    tri_facets.append_unique( facet3 );
    tri_facets.append_unique( facet4 );
  }
  else if( facet1 )
  {
    tri_facets.append_unique( facet1 );
    tri_facets.append( new CubitFacetData( pt3, pt4, pt1 ) );
  }
  else if( facet2 )
  {
    tri_facets.append_unique( facet2 );
    tri_facets.append( new CubitFacetData( pt1, pt2, pt3 ) );
  }
  else if( facet3 )
  {
    tri_facets.append_unique( facet3 );
    tri_facets.append( new CubitFacetData( pt2, pt3, pt4 ) );
  }
  else if( facet4 )
  {
    tri_facets.append_unique( facet4 );
    tri_facets.append( new CubitFacetData( pt1, pt2, pt4 ) );
  }
  else
  {
    tri_facets.append( new CubitFacetData( pt1, pt2, pt3 ) );
    tri_facets.append( new CubitFacetData( pt3, pt4, pt1 ) );
  }
}


void ExoLoader::eval_surface( Surface* surf_ptr )
{
  RefFace* face = dynamic_cast<RefFace*>(surf_ptr->topology_entity());
  FacetSurface* facet_surf_ptr = dynamic_cast<FacetSurface*> (surf_ptr);
  
  if( face )
    PRINT_INFO("RefFace %d:\n",face->id());
  else 
    PRINT_INFO("Surface @%p:\n",surf_ptr);
    
  const char* type;
  switch( surf_ptr->geometry_type() )
  {
    case CONE_SURFACE_TYPE:
      type = "conic";
      break;
    case PLANE_SURFACE_TYPE:
      type = "planar";
      break;
    case SPHERE_SURFACE_TYPE:
      type = "spherical";
      break;
    case SPLINE_SURFACE_TYPE:
      type = "spline";
      break;
    case TORUS_SURFACE_TYPE:
      type = "toroidal";
      break;
    case BEST_FIT_SURFACE_TYPE:
      type = "best fit";
      break;
    case FACET_SURFACE_TYPE:
      type = "facet";
      break;
    default:
      type = "unknown";
  }
  PRINT_INFO("\tSurface type: %s\n",type);
  
  PRINT_INFO("\tArea:         %lf\n",surf_ptr->measure());
  
  CubitVector p1, p2, p3;
  
  if( surf_ptr->get_point_normal( p1, p2 ) )
  {
    PRINT_INFO("\tPlanar:       YES\n"
               "\tNormal:       (%lf,%lf,%lf) at (%lf,%lf,%lf)\n",
               p2.x(),p2.y(),p2.z(),p1.x(),p1.y(),p1.z());
  }
  else
  {
    PRINT_INFO("\tPlanar:       NO\n");
  }
  
  CubitBox bbox = surf_ptr->bounding_box();
  PRINT_INFO("\tBounding Box: X:[% lf, % lf],\n"
             "\t              Y:[% lf, % lf],\n"
             "\t              Z:[% lf, % lf]\n",
             bbox.minimum().x(), bbox.maximum().x(),
             bbox.minimum().y(), bbox.maximum().y(),
             bbox.minimum().z(), bbox.maximum().z() );
             
  DLIList<CubitFacet*> facet_list;
  DLIList<CubitPoint*> point_list;
  if(facet_surf_ptr){
    
    if( facet_surf_ptr->get_my_facets( facet_list, point_list ) )
    {
      PRINT_INFO("\tFacet count:   %d\n", facet_list.size() );
      PRINT_INFO("\tPoint count:   %d\n", point_list.size() );
    }
    else
    {
      PRINT_ERROR("FacetSurface::get_my_facets(..) returned CUBIT_FAILURE.\n");
    } 
  }
  else
    PRINT_ERROR(" Surface was not a facet surface. \n");
  
  
  CubitVector points[9];
  bbox.get_corners( points );
  points[8] = bbox.center();
  
  PRINT_INFO("\t  Source Point    Closest Point   Normal Closest  Closest Trimmed\n"
             "\t  --------------  --------------  --------------  ---------------\n");
  for( int i = 0; i < 9; i++ )
  {
    CubitStatus s = surf_ptr->closest_point( points[i], &p1, &p2 );
    surf_ptr->closest_point_trimmed( points[i], p3 );
    
    if( s )
    {
      PRINT_INFO("\tX:% 14.6lf  % 14.6lf  % 14.6lf  % 14.6lf\n",
                   points[i].x(), p1.x(),   p2.x(),   p3.x() );
      PRINT_INFO("\tY:% 14.6lf  % 14.6lf  % 14.6lf  % 14.6lf\n",
                   points[i].y(), p1.y(),   p2.y(),   p3.y() );
      PRINT_INFO("\tZ:% 14.6lf  % 14.6lf  % 14.6lf  % 14.6lf\n",
                   points[i].z(), p1.z(),   p2.z(),   p3.z() );
    }
    else
    {
      PRINT_INFO("\tX:% 13.6lf  FAILED          FAILED          FAILED\n", 
        points[i].x() );
      PRINT_INFO("\tY:% 13.6lf  FAILED          FAILED          FAILED\n", 
        points[i].y() );
      PRINT_INFO("\tZ:% 13.6lf  FAILED          FAILED          FAILED\n",
        points[i].z() );
    }
    PRINT_INFO("\n");
  }
  
  
  
  if( ! surf_ptr->is_parametric() )
  {
    PRINT_INFO("\tParametric:   NO.\n");
    return;
  }
  
  double u1, u2, v1, v2;
  if( ! surf_ptr->get_param_range_U( u1, u2 ) )
  {
    PRINT_ERROR("Surface claims to be parametric, but "
                "get_param_range_U(..) returned failuer.\n");
  }
  if( ! surf_ptr->get_param_range_V( v1, v2 ) )
  {
    PRINT_ERROR("Surface claims to be parametric, but "
                "get_param_range_V(..) returned failuer.\n");
  }
  PRINT_INFO("\tParam Range:  U:[%lf,%lf], V:[%lf,%lf]\n", u1, u2, v1, v2 );
  
  if( surf_ptr->is_periodic() )
  {
    double p;
    if( surf_ptr->is_periodic_in_U(p) )
      PRINT_INFO("\tU-Period:     %lf\n", p );
    if( surf_ptr->is_periodic_in_V(p) )
      PRINT_INFO("\tV-Period:     %lf\n", p );
  }
  else
  {
    PRINT_INFO("\tPeriodic:     NO\n");
  }
  
  if( surf_ptr->is_closed_in_U() )
    PRINT_INFO("\tSurface is closed in U direction.\n");
  if( surf_ptr->is_closed_in_V() )
    PRINT_INFO("\tSurface is closed in V direction.\n");
  
}
