
#include "FacetEntity.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "CubitFacetEdge.hpp"
#include "GfxDebug.hpp"
#include "CubitBox.hpp"

#include "debug.hpp"
static int fg_color = CUBIT_MAGENTA;
void dcolor(int icol)
{
//const int CUBIT_DEFAULT_COLOR = -1;
//const int CUBIT_BLACK         =  0;
//const int CUBIT_GREY          =  1;
//const int CUBIT_GREEN         =  2;
//const int CUBIT_YELLOW        =  3;
//const int CUBIT_RED           =  4;
//const int CUBIT_MAGENTA       =  5;
//const int CUBIT_CYAN          =  6;
//const int CUBIT_BLUE          =  7;
//const int CUBIT_WHITE         =  8;
//const int CUBIT_ORANGE        =  9;
//const int CUBIT_PINK          = 10;
//const int CUBIT_BROWN         = 11;
//const int CUBIT_GOLD          = 12;
//const int CUBIT_LIGHTBLUE     = 13;
//const int CUBIT_LIGHTGREEN    = 14;
//const int CUBIT_SALMON        = 15;
  fg_color = icol;
}

void ddraw( FacetEntity *facet_ptr )
{
  facet_ptr->debug_draw( fg_color );
}

void dfdraw( CubitFacet *facet_ptr  )
{
  facet_ptr->debug_draw( fg_color );
}

void dedraw( CubitFacetEdge *facet_ptr )
{
  facet_ptr->debug_draw( fg_color );
}

void dpdraw( CubitPoint *facet_ptr )
{
  facet_ptr->debug_draw( fg_color );
}


void dview()
{
  GfxDebug::flush();
  GfxDebug::mouse_xforms();
}

void dzoom(CubitBox &box)
{
  GfxDebug::zoom(box);
}

void dldraw( DLIList<FacetEntity *> &facet_list )
{
  FacetEntity *facet_ptr;
  for (int ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->debug_draw( fg_color, 0 );
  }
  GfxDebug::flush();
}

void dfldraw( DLIList<CubitFacet *> &facet_list )
{
  CubitFacet *facet_ptr;
  for (int ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->debug_draw( fg_color, 0 );
  }
  GfxDebug::flush();
}

void deldraw( DLIList<CubitFacetEdge *> &edge_list )
{
  CubitFacetEdge *edge_ptr;
  for (int ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    edge_ptr->debug_draw( fg_color, 0 );
  }
  GfxDebug::flush();
}

void dpldraw( DLIList<CubitPoint *> &point_list )
{
  CubitPoint *point_ptr;
  for (int ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    point_ptr->debug_draw( fg_color, 0 );
  }
  GfxDebug::flush();
}

int dflcheck( DLIList<CubitFacet *> &facet_list )
{
  int ii;
  int ier = 0;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    CubitFacet *facet_ptr = facet_list.get_and_step();
    ier += dfcheck( facet_ptr );
  }
  return ier;
}

int dcheck( DLIList<FacetEntity *> &facet_list )
{
  int ii;
  int ier = 0;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    CubitFacet *facet_ptr = (CubitFacet *)facet_list.get_and_step();
    ier += dfcheck( facet_ptr );
  }
  return ier;
}

int dfcheck( CubitFacet *facet_ptr )
{
  int ier = 0;
  
  // check edges

  int ii, jj;
  for (ii=0; ii<3; ii++)
  {
    CubitFacetEdge *edge_ptr = facet_ptr->edge( ii );
    DLIList <CubitFacet *> facet_list;
    edge_ptr->facets( facet_list );
    int found = 0;
    for (jj=0; jj<facet_list.size() && !found; jj++)
    {
      CubitFacet *f = facet_list.get_and_step();
      if (f == facet_ptr)
      {
        found = 1;
      }
    }
    if (!found)
    {
      PRINT_ERROR( "Facet %d is not in Edge %d adjacency list\n",
                    facet_ptr->id(), edge_ptr->id() );
      ier++;
    }

    // check the edge's orientation on the facet

    CubitPoint *ep0, *ep1;
    int use = facet_ptr->edge_use( ii );
    if (use > 0)
    {
      ep0 = edge_ptr->point( 0 );
      ep1 = edge_ptr->point( 1 );
    }
    else
    {
      ep1 = edge_ptr->point( 0 );
      ep0 = edge_ptr->point( 1 );
    }

    CubitPoint *fp0 = facet_ptr->point( (ii+1)%3 );
    CubitPoint *fp1 = facet_ptr->point( (ii+2)%3 );
    if (fp0 != ep0 || fp1 != ep1)
    {
      PRINT_ERROR( "Edge %d on Facet %d is not oriented with points %d and %d correctly\n",
        edge_ptr->id(), facet_ptr->id(), fp0->id(), fp1->id() );
      ier++;
    }
  }

  // check nodes 

  for (ii=0; ii<3; ii++)
  {
    CubitPoint *point_ptr = facet_ptr->point( ii );
    DLIList <CubitFacet *> facet_list;
    point_ptr->facets( facet_list );
    int found = 0;
    for (jj=0; jj<facet_list.size() && !found; jj++)
    {
      CubitFacet *f = facet_list.get_and_step();
      if (f == facet_ptr)
      {
        found = 1;
      }
    }
    if (!found)
    {
      PRINT_ERROR( "Facet %d is not in Point %d adjacency list\n",
                    facet_ptr->id(), point_ptr->id() );
      ier++;
    }
  }
  return ier;
}

// draw a vector
void dray( const CubitVector &start, const CubitVector &vec, double length )
{
  CubitVector end = start+length*vec;
  GfxDebug::draw_vector(end, start, fg_color );
  GfxDebug::flush();
}

// draw a point
void dpoint( const CubitVector &pt )
{
  GfxDebug::draw_point(pt, fg_color);
  GfxDebug::flush();

}

int get_color()
{
  return fg_color;
}

//EOF

