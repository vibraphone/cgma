#include "PST_Data.hpp"
#include "GMem.hpp"
#include "GfxDebug.hpp"
#include "GeometryDefines.h"
#define PST_MAX_LIST_LEN 1000
// Used only for debugging.  Has no effect
// if compiled -DNDEBUG.

const double RESABS_SQR = GEOMETRY_RESABS * GEOMETRY_RESABS;


/*************************************************************************
 *                               PST_Point
 ************************************************************************/

PST_Point::~PST_Point()
{
  while( edge_ )
    delete edge_;
}

PST_Edge* PST_Point::common( PST_Point* pt )
{
  if( edge_ )
  {
    int debug = 0;
    
    PST_Edge* e = edge_;
    do
    {
      if( e->other(pt) == this )
        return e;
      e = next( e );  
    
      assert( ++debug < PST_MAX_LIST_LEN );
    } while( e && e != edge_ );
  }
  return 0;
}


int PST_Point::faces( DLIList<PST_Face*>* list )
{
  int count = 0;
  if( edge() )
  {
    PST_Edge* e = edge();
    do
    {
      if (e->forward()->face())
        e->forward()->face()->private_mark_ = true;
      if (e->reverse()->face())
        e->reverse()->face()->private_mark_ = true;
      e = e->next(this);
    } while( e != edge() );
    
    do
    {
      if( e->forward()->face() && e->forward()->face()->private_mark_ )
      {
        count++;
        e->forward()->face()->private_mark_ = false;
        if( list ) list->append(e->forward()->face());
      }
      if( e->reverse()->face() && e->reverse()->face()->private_mark_ )
      {
        count++;
        e->reverse()->face()->private_mark_ = false;
        if( list ) list->append(e->reverse()->face());
      }
      e = e->next(this);
    } while( e != edge() );
  }      
  
  return count;
}

/*************************************************************************
 *                               PST_Face
 ************************************************************************/

void PST_Face::append_points( DLIList<PST_Point*>& result_list )
{
  PST_CoEdge* ce = coedge_;
  int debug = 0;
  do
  {
    result_list.append( ce->end_point() );
    ce = ce->next();
    
    assert( ++debug < PST_MAX_LIST_LEN );
  } while( ce != coedge_ );
}

PST_Point* PST_Face::opposite( PST_Edge* edge )
{
  PST_CoEdge* ce = coedge_;
  int debug = 0;
  do
  {
    if( ce->edge() == edge )
      return ce->next()->end_point();
    ce = ce->next();
    
    assert( ++debug < PST_MAX_LIST_LEN );
  } while( ce != coedge_ );
  
  return 0;
}

PST_Edge* PST_Face::opposite( PST_Point* point )
{
  PST_CoEdge* ce = coedge_;
  int debug = 0;
  do
  {
    if( ce->start_point() == point )
      return ce->next()->edge();
    ce = ce->next();
    
    assert( ++debug < PST_MAX_LIST_LEN );
  } while( ce != coedge_ );
  
  return 0;
}

bool PST_Face::two_edges( PST_Point* point, PST_Edge*& e1, PST_Edge*& e2 )
{
  e1 = e2 = 0;
  PST_CoEdge* ce = coedge_;
  int debug = 0;
  do
  {
    if( ce->end_point() == point )
    {
      if( e1 || e2 )
        return false;
      e1 = ce->edge();
      e2 = ce->next()->edge();
    }
    ce = ce->next();
    
    assert( ++debug < PST_MAX_LIST_LEN );
  } while( ce != coedge_ );
  
  return e1 && e2;
}


bool PST_Face::calculate_plane()
{
  update_plane_ = 0;

  // Use Newell's method
  
  CubitVector dif, sum, ref(0.,0.,0.);
  CubitVector norm(0.,0.,0.);
  
  // For each coedge...
  
  PST_CoEdge* ce = coedge_;
  int count = 0;
  do
  {
    count++;
    const CubitVector& pt1 = ce->start_coord();
    const CubitVector& pt2 = ce->end_coord();
    dif   = pt2 - pt1;
    sum   = pt2 + pt1;
    ref  += pt1;
    norm += CubitVector( dif.y()*sum.z(), dif.z()*sum.x(), dif.x()*sum.y() );
    ce = ce->next();
    
    assert( count < PST_MAX_LIST_LEN );
  } while( ce != coedge_ );
  
  // Degenerate?
  
  double len = norm.length();
  if( len > CUBIT_RESABS )
  {
    double d = (ref % norm) / (count * len);
    norm /= -len; //reverse and normalize

    plane_.normal( norm );
    plane_.coefficient( d );

    return true;
  }
  else
  {
    plane_.normal( CubitVector(0.,0.,0.) );
    plane_.coefficient( 0. );

    return false;
  }
}      
    


/*************************************************************************
 *                              PST_CoEdge
 ************************************************************************/
PST_CoEdge* PST_CoEdge::previous()
{
  if( !next_ )
    return 0;
    
  PST_CoEdge* coedge = this;
  int debug = 0;
  while( coedge->next_ != this )
  {
    coedge = coedge->next_;
    assert( ++debug < PST_MAX_LIST_LEN );
  }
  return coedge;
}

/*************************************************************************
 *                               PST_Edge
 ************************************************************************/

PST_Edge::~PST_Edge()
{
    // Update list "head" pointers on points, if necessary.
  if( start_point() )
    remove_start_point();
  if( end_point() )
    remove_end_point();

    // Destroy parent facets
  for ( int i = 0; i < 2; i++ )
  {
    PST_CoEdge* coedge = i ? reverse() : forward();
    if ( coedge->face() )
    {
      PST_Face* face = coedge->face();
      face->coedge_ = 0;
      
      PST_CoEdge* first = coedge;
      do {
        assert(coedge->face() == face);
        coedge->face_ = 0;
        PST_CoEdge* next = coedge->next();
        coedge->next_ = 0;
        coedge = next;
      } while( first != coedge );
      delete face ;
    }
  }
}
    
double PST_Edge::closest_on_line( const CubitVector& P )
{
  CubitVector B = start_coord();
  CubitVector M   = end_coord() - B;
  if( M.length_squared() < RESABS_SQR )
    return 0.0;
  return ( M % ( P - B ) ) / ( M % M );
}

double PST_Edge::closest_on_edge( const CubitVector& p )
{
  double t = closest_on_line( p );
  if( t < 0.0 )
    t = 0.0;
  else if( t > 1.0 )
    t = 1.0;
  return t;
}

double PST_Edge::closest_on_line( const CubitVector& B2, 
                                  const CubitVector& M2 )
{
  CubitVector B1 = start_coord();
  CubitVector M1 = direction();
  
  if( M1.length_squared() < RESABS_SQR )
    return 0.0;
  
  if( M2.length_squared() < RESABS_SQR )
    return closest_on_line( B2 );
  
  CubitVector cross = M2 * M1;
  if( cross.length_squared() < CUBIT_RESABS ) //parallel
    return 0.0;
  
  CubitVector N = M2 * cross;
  double      D = -( N % B2 );
  return -( N % B1 + D ) / ( N % M1 );
}


PST_Edge* PST_Edge::other( PST_Point* point, PST_Face* face )
{
  PST_Edge *result = 0;
  for( PST_Edge* e = point->next(this); e != this; e = point->next(e) )
  {
    if( e->forward()->face() == face || e->reverse()->face() == face )
    {
      if( result && result != e )
        return 0;
      result = e;
    }
  }
  return result;
}

/*************************************************************************
 *                         Construction Methods
 ************************************************************************/

PST_Face* PST_Edge::create_face( PST_Point* pt1, PST_Point* pt2, PST_Point* pt3 )
{
    // It is trivial to generalize this function for creating
    // faces with an arbitrary number of sides, but I just
    // want triangles.
  const int size = 3;
  PST_Point* p[size] = { pt1, pt2, pt3 };

  int i;
  PST_Edge* e[size];
  PST_CoEdge* c[size];
  
    // Initialize edge and coedge arrays, creating
    // edges if necessary.
  for( i = 0; i < size; i++ )
  {
    e[i] = p[i]->common( p[(i+1)%size] );
    if( !e[i] )
    {
      e[i] = new PST_Edge( p[i], p[(i+1)%size] );
      c[i] = e[i]->forward();
    }
    else
    {
      c[i] = e[i]->start_point() == p[i] ? e[i]->forward() : e[i]->reverse();
      
        // If the co-edge in the appropriate direction is already
        // part of a face (other than the boundary face), then we
        // can't proceed.
//      if( !c[i]->face()->boundary() )
      if( c[i]->face() )
      {
        while( --i >= 0 )
          if( ! c[i]->face() )
            delete e[i];
        return 0;
      }
    }
  }
  
  PST_Face* result = new PST_Face( c[0] );
  for ( i = 0; i < size; i++ )
  {
    int next = (i + 1) % size;
    c[i]->face_ = result;
    c[i]->next_ = c[next];
  }
  return result;

}


PST_Edge* PST_Edge::split_face( PST_Point* start, PST_Point* end, PST_Face* face )
{
  if( !start->edge_ || !end->edge_ )
    return 0;
  
  if( ! face )
  {
    PST_Edge* start_edge = start->edge();
    do
    {
      PST_Face* curr_face = start_edge->forward()->face();
      PST_Face* rev_face = start_edge->reverse()->face();
      while( true )
      {
        PST_CoEdge* coe = curr_face->coedge_;
        do
        {
          if( coe->other( end ) )
          {
            if( face ) return 0;
            
            face = curr_face;
            curr_face = rev_face;
            break;
          }
          coe = coe->next();
        } while( coe != curr_face->coedge_ );
        
        if( curr_face == rev_face )
          break;
        else
          curr_face = rev_face;
      }

      start_edge = start_edge->next(start);
    } while( start_edge != start->edge() );
    
    if( !face ) return 0;
  }
  
  PST_CoEdge* start_coedge = 0;
  PST_CoEdge*   end_coedge = 0;
  
  PST_CoEdge* coedge = face->coedge_;
  do
  {
    if( coedge->end_point() == start )
    {
      start_coedge = coedge;
      break;
    }
    coedge = coedge->next();
  } while( coedge != face->coedge_ );
  
  if( ! start_coedge ) 
    return 0;
  
  coedge = start_coedge->next();
  PST_CoEdge* stop = start_coedge;
  do
  {
    if( coedge->end_point() == start )
    {
      start_coedge = coedge;
    }
   
    if( coedge->end_point() == end )
    {
      end_coedge = coedge;
      break;
    }
    
    coedge = coedge->next();
  } while( coedge != stop );
  
  if( ! end_coedge )
    return 0;
  
  
  return split_face( start_coedge, end_coedge );
}

PST_Edge* PST_Edge::split_face( PST_CoEdge* coedge1, PST_CoEdge* coedge2 )
{
  PST_Point* point1 = coedge1->end_point();
  PST_Point* point2 = coedge2->end_point();
  
  if( (coedge1->edge()      == coedge2->edge() ) ||
      (coedge1->face()      != coedge2->face() ) ||
      (point1               == point2          ) ||
      ( point1->common( point2 )               ) )
    return 0;
  
  PST_Edge* new_edge = new PST_Edge( point1, point2 );
  new_edge->forward_.next_ = coedge2->next_;
  new_edge->reverse_.next_ = coedge1->next_;
  coedge1->next_ = new_edge->forward();
  coedge2->next_ = new_edge->reverse();
  
  PST_Face* old_face = coedge1->face();
  PST_Face* new_face = new PST_Face( new_edge->forward(), old_face );
  //new_face->boundary(0);
  old_face->coedge_ = new_edge->reverse();
  new_edge->reverse_.face_ = old_face;
  
  coedge1 = new_edge->forward();
  do
  {
    coedge1->face_ = new_face;
    coedge1 = coedge1->next();
  } while( coedge1 != new_edge->forward() );

  return new_edge;
}


PST_Edge* PST_Edge::insert_in_face( PST_Point* end, 
                                    PST_CoEdge* after_this )
{
  if( end->edge_ )
    return 0;
  
  PST_Point* start = after_this->end_point();
  PST_Edge* new_edge = new PST_Edge( start, end );
  new_edge->forward_.next_ = new_edge->reverse();
  new_edge->reverse_.next_ = after_this->next_;
  new_edge->forward_.face_ = new_edge->reverse_.face_ = after_this->face_;
  after_this->next_ = new_edge->forward();
  after_this->face_->modified(true);
  return new_edge;
}


PST_Edge* PST_Edge::split( PST_Point* point )
{
  if( point->edge_ )
    return 0; 
  
  PST_Edge* new_edge = new PST_Edge( point, end_ );

  if (forward()->face()) {
    new_edge->forward_.face_ = forward_.face_;
    new_edge->forward_.next_ = forward_.next_;
    forward_.next_ = new_edge->forward();
    forward_.face_->modified(true);
  }
  if (reverse()->face()) {
    new_edge->reverse_.face_ = reverse_.face_;
    new_edge->reverse_.next_ = reverse();
    reverse_.previous()->next_ = new_edge->reverse();
    reverse_.face_->modified(true);
  }  

  set_end_point(point);
  
  return new_edge;
}


void PST_Edge::set_point( PST_Point* pt, 
                          PST_Point*& ptr,
                          PST_Edge*& next )
{
  if( ptr )
    remove_point( ptr, next );
  
  
  if( pt->edge_ )
  {
    if( pt->edge()->start_point() == pt )
    {
      next = pt->edge()->start_next_;
      pt->edge()->start_next_ = this;
    }
    else
    {
      assert( pt->edge()->end_point() == pt );
      next = pt->edge()->end_next_;
      pt->edge()->end_next_ = this;
    }
  }
  else
  {
    next = this;
    pt->edge_ = this;
  }
    
  ptr = pt;
}

void PST_Edge::remove_point( PST_Point*& ptr, PST_Edge*& next )
{
  if( next == this )
  {
    ptr->edge_ = 0;
    next= 0;
    ptr = 0;
    return;
  }
  
  PST_Edge* prev = next;
  while( prev->next(ptr) != this )
  {
    prev = prev->next(ptr);
    assert( prev != next );
  }
  
  if( prev->start_next_ == this )
  {
    prev->start_next_ = next;
  }
  else
  {
    assert(prev->end_next_ == this);
    prev->end_next_ = next;
  }
  
  if( ptr->edge_ == this )
    ptr->edge_ = next;
}


void PST_Point::debug_draw( int color, bool flush )
{
  GfxDebug::draw_point( float(x()), float(y()), float(z()), color );
  if( flush ) 
    GfxDebug::flush();
}
  
void PST_Edge::debug_draw( int color, bool flush )
{
  GfxDebug::draw_line( 
    float(start_coord().x()), float(start_coord().y()), float(start_coord().z()),
    float(  end_coord().x()), float(  end_coord().y()), float(  end_coord().z()),
    color );
  if( flush ) GfxDebug::flush();
}

void PST_Face::debug_draw( int color, bool flush )
{
  PST_CoEdge* coedge = first_coedge();
  do
  {
    coedge->edge()->debug_draw( color, false );
    coedge = coedge->next();
  } while( coedge != first_coedge() );
  if( flush ) GfxDebug::flush();
}

void PST_Edge::make_gmem( GMem& gmem, DLIList<PST_Face*>& facets )
{
  DLIList<PST_Point*> point_list, temp_list;
  int i, j, pcount, fcount;
  
  for( i = facets.size(); i--; )
  {
    PST_Face* facet = facets.get_and_step();
//    if( facet->boundary() )
//      continue;

    temp_list.clean_out();
    facet->append_points( temp_list );
    for( j = temp_list.size(); j--; )
    {
      temp_list.get_and_step()->mark = -1;
    }
  }
  
  pcount = 0;
  fcount = 0;
  for( i = facets.size(); i--; )
  {
    PST_Face* facet = facets.get_and_step();
//    if( facet->boundary() )
//      continue;
    
    temp_list.clean_out();
    facet->append_points( temp_list );
    fcount += temp_list.size() + 1;
    for( j = temp_list.size(); j--; )
    {
      PST_Point* pt_ptr = temp_list.get_and_step();
      if( pt_ptr->mark == -1 )
      {
        pt_ptr->mark = pcount++;
        point_list.append( pt_ptr );
      }
    }
  }
  
  gmem.allocate_tri( facets.size() );
  gmem.pointListCount = point_list.size();
  GPoint* pt_array = gmem.point_list();
  
  point_list.reset();
  for( i = 0; i < point_list.size(); i++ )
  {
    PST_Point* pt_ptr = point_list.get_and_step();
    pt_array[i].x = float(pt_ptr->x());
    pt_array[i].y = float(pt_ptr->y());
    pt_array[i].z = float(pt_ptr->z());
  }
  
  gmem.fListCount = facets.size() * 4;
  int* offset = gmem.facet_list();
  for( i = facets.size(); i--; )
  {
//    if( facets.get()->boundary() )
//      continue;
      
    temp_list.clean_out();
    facets.get_and_step()->append_points( temp_list );
    *(offset++) = temp_list.size();
    for( j = temp_list.size(); j--; )
    {
      *(offset++) = temp_list.get_and_step()->mark;
    }
  }
}

void PST_Edge::make_gmem( GMem& gmem, DLIList<PST_Edge*>& edges )
{
  DLIList<PST_Face*> face_list;
  PST_Edge::faces( edges, face_list );
  make_gmem( gmem, face_list );
}

void PST_Point::debug_draw_points( DLIList<PST_Point*>& point_list, int color, bool flush )
{
  for( int p = point_list.size(); p--; )
    point_list.get_and_step()->debug_draw(color, false);
  if( flush ) GfxDebug::flush();
}

void PST_Edge::debug_draw_points( DLIList<PST_Edge*>& edge_list, 
                            int color, int boundary_color, bool flush )
{
  if( ! boundary_color ) boundary_color = color;
  
  DLIList<PST_Point*> point_list;
  DLIList<PST_Point*> boundary_list;
  PST_Edge::points( edge_list, point_list );
  for( int i = point_list.size(); i--; )
  {
     PST_Point* pt = point_list.step_and_get();
     PST_Edge* e = pt->edge();
     do
     {
       if( !e->forward()->face() || !e->reverse()->face() )
         break;
       e = e->next(pt);
     } while( e != pt->edge() );
     if( !e->forward()->face() || !e->reverse()->face() )
     {
       point_list.change_to( 0 );
       boundary_list.append( pt );
     }
  }
  point_list.remove_all_with_value(0);
  
  PST_Point::debug_draw_points( point_list, color, false );
  PST_Point::debug_draw_points( boundary_list, boundary_color, flush );
}

void PST_Edge::debug_draw_edges(  DLIList<PST_Edge*>& edge_list,
                            int color, int boundary_color,
                            bool flush )
{
  if( !boundary_color ) boundary_color = color;
  
  for( int e = edge_list.size(); e--; )
  {
    PST_Edge* edge_ptr = edge_list.get_and_step();
    int c = !edge_ptr->forward()->face() || !edge_ptr->reverse()->face()
       ? boundary_color : color;
    edge_ptr->debug_draw( c, false );
  }
  if( flush ) GfxDebug::flush();
}

void PST_Edge::debug_draw_faces( DLIList<PST_Edge*>& edge_list, int color, bool flush )
{
  DLIList<PST_Face*> face_list;
  PST_Edge::faces( edge_list, face_list );
  PST_Face::debug_draw_faces( face_list, color, flush );
}

void PST_Edge::faces( DLIList<PST_Edge*>& edges, DLIList<PST_Face*>& faces )
{
  int e;
  for( e = edges.size(); e--; )
  {
    PST_Edge* edge_ptr = edges.get_and_step();
    if( edge_ptr->forward()->face() )
      edge_ptr->forward()->face()->private_mark_ = 1;
    if( edge_ptr->reverse()->face() )
      edge_ptr->reverse()->face()->private_mark_ = 1;
  }
  for( e = edges.size(); e--; )
  {
    PST_Edge* edge_ptr = edges.get_and_step();
    PST_Face* fface_ptr = edge_ptr->forward()->face();
    PST_Face* rface_ptr = edge_ptr->reverse()->face();
    if( fface_ptr && fface_ptr->private_mark_ )
    {
      fface_ptr->private_mark_ = 0;
      faces.append( fface_ptr );
    }
    if( rface_ptr && rface_ptr->private_mark_ )
    {
      rface_ptr->private_mark_ = 0;
      faces.append( rface_ptr );
    }
  }
}

void PST_Edge::edges( DLIList<PST_Face*>& faces, DLIList<PST_Edge*>& edges )
{
  int f;
  for( f = faces.size(); f--; )
  {
    PST_CoEdge* first = faces.get_and_step()->first_coedge();
    PST_CoEdge* coedge = first;
    do
    {
      coedge->edge()->private_mark_ = 1;
      coedge = coedge->next();
    } while( coedge != first );
  }
  
  for( f = faces.size(); f--; )
  {
    PST_CoEdge* first = faces.get_and_step()->first_coedge();
    PST_CoEdge* coedge = first;
    do
    {
      if( coedge->edge()->private_mark_ )
      {
        coedge->edge()->private_mark_ = 0;
        edges.append( coedge->edge() );
      }
      coedge = coedge->next();
    } while( coedge != first );
  }
}  

void PST_Edge::edges( DLIList<PST_Point*>& pts, DLIList<PST_Edge*>& edges )
{
  int p;
  for( p = pts.size(); p--; )
  {
    PST_Point* pt = pts.get_and_step();
    PST_Edge* edge = pt->edge();
    if( edge ) do
    {
      edge->private_mark_ = 1;
      edge = edge->next( pt );
    } while( edge != pt->edge() );
  }
  
  for( p = pts.size(); p--; )
  {
    PST_Point* pt = pts.get_and_step();
    PST_Edge* edge = pt->edge();
    if( edge ) do
    {
      if( edge->private_mark_ )
      {
        edge->private_mark_ = 0;
        edges.append( edge );
      }
      edge = edge->next( pt );
    } while( edge != pt->edge() );
  }
}


void PST_Edge::points( DLIList<PST_Edge*>& edges, DLIList<PST_Point*>& points )
{
  int e;
  for( e = edges.size(); e--; )
  {
    PST_Edge* edge_ptr = edges.get_and_step();
    edge_ptr->start_point()->private_mark_ = 1;
    edge_ptr->end_point()->private_mark_ = 1;
  }
  for( e = edges.size(); e--; )
  {
    PST_Edge* edge_ptr = edges.get_and_step();
    PST_Point* sp = edge_ptr->start_point();
    PST_Point* ep = edge_ptr->end_point();
    if( sp->private_mark_ )
    {
      sp->private_mark_ = 0;
      points.append( sp );
    }
    if( ep->private_mark_ )
    {
      ep->private_mark_ = 0;
      points.append( ep );
    }
  }
}

void PST_Face::debug_draw_faces( DLIList<PST_Face*>& face_list, int color, bool flush )
{
  for( int f = face_list.size(); f--; )
    face_list.get_and_step()->debug_draw( color, false );
  if( flush ) GfxDebug::flush();
}

void PST_Edge::make_facets( GMem& gmem, double tolerance, 
                            DLIList<PST_Edge*>& edge_list )
{
  assert(gmem.fListCount % 4 == 0);
  std::vector<double> points(gmem.pointListCount*3);
  std::vector<int> facets(gmem.fListCount*3/4);
  int i;
  GPoint* pitor = gmem.point_list();
  std::vector<double>::iterator ditor = points.begin();
  for ( i = gmem.pointListCount; i--; )
  {
    *ditor++ = pitor->x;
    *ditor++ = pitor->y;
    *ditor++ = pitor->z;
    pitor++;
  }
  
  int* fitor = gmem.facet_list();
  std::vector<int>::iterator iitor = facets.begin();
  for ( i = 0; i < gmem.fListCount; i += 4 )
  {
    assert( *fitor++ == 3 );
    *iitor++ = *fitor++;
    *iitor++ = *fitor++;
    *iitor++ = *fitor++;
  }
   
  make_facets( points, facets, tolerance, edge_list ); 
}

void PST_Edge::make_facets( 
    const std::vector<double>& coordinates,
    const std::vector<int>& connections,
    double tolerance, 
    DLIList<PST_Edge*>& edge_list )
{
  DLIList<PST_Face*> face_list;
  double tol_sqr = tolerance * tolerance;
  int i, j, k, numcoords, numtriangles;
  
  numcoords = coordinates.size()/3;
  numtriangles = connections.size()/3;
  //The list of points created.
  PST_Point** ptlist = new PST_Point*[numcoords];
  
  //The list of indices into ptlist, where the index into this
  //list is the same as the index into the point list, and
  //the value in this list is the index of the corresponding 
  //point in ptlist.  This is used because some points will
  //be merged due to the tolerance.  
  int* ptindex_list = new int[numcoords];
   
  //Create the points
  k = 0;
  for( i = 0; i < numcoords; i++ )
  {
    CubitVector pos( coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2] );
 
    //Check if we need to merge this point with a different point
    for( j = 0; j < k; j++ )
    {
      if( (*(ptlist[j]) - pos).length_squared() <= tol_sqr )
        break;
    }
    if( j == k )
    {
      ptlist[k] = new PST_Point( pos );
      ptlist[k++]->sequence = i;
    }
    ptindex_list[i] = j;
  } 
  
  int fail_count = 0;
  for ( i = 0; i < numtriangles; i++ )
  {
    int i1 = ptindex_list[connections[3*i  ]];
    int i2 = ptindex_list[connections[3*i+1]];
    int i3 = ptindex_list[connections[3*i+2]];
    if ( i1 == i2 || i2 == i3 || i3 == i1 )
    {
      PRINT_ERROR("Degenerate facet encountered in PST_Edge::make_facets()\n");
      fail_count++;
    }
    else
    {
      PST_Face* new_face 
        = PST_Edge::create_face( ptlist[i1], ptlist[i2], ptlist[i3] );
      if( new_face ){
        face_list.append(new_face);
        new_face->sequence = i;
      }
      else{
        fail_count++;
      }
    }
  }
   
  if( fail_count > 0 )
    PRINT_ERROR("Failed to construct %d facets in PST_Data::make_facets(..)\n",
      fail_count);
   
  delete [] ptindex_list;
  delete [] ptlist;
  
  if(fail_count == 0)
  {
    PST_Edge::edges( face_list, edge_list );
    validate( edge_list, true );
    bool debug1 = false;
    if (debug1)
      debug_draw_edges( edge_list, CUBIT_BLUE, CUBIT_RED, true );
  }
}

int PST_Point::validate(CubitBoolean print)
{
  if( edge_ && !edge_->other(this) )
  {
    if( print )
    {
      PRINT_INFO("Bad point->edge link.  Point %p, Edge %p.\n",
                 static_cast<void*>(this), static_cast<void*>(edge_));
    }
    return 1;
  }
  return 0;
}

int PST_Edge::validate(CubitBoolean print)
{
  int result = 0;
  int count = 0;
  PST_Edge* edge = 0;

  if( !start_ )
  {
    if( print ) PRINT_ERROR("Edge %p has null start point.\n",
                            static_cast<void*>(this));
    result++;
  }
  else 
  {
    result += start_->validate(print);
  
    edge = this;
    count = 0;
    do
    {
      if( count++ > PST_MAX_LIST_LEN )
      {
        if( print )
        {
          PRINT_ERROR("Bad edge list around Point %p.\n",
                      static_cast<void*>(start_));
        }
        result++;
        break;
      }
      edge = edge->next( start_ );
    } while( edge != this );
  }  

  
  if( !end_ )
  {
    if( print ) PRINT_ERROR("Edge %p has null end point.\n",
                            static_cast<void*>(this));
    result++;
  }
  else 
  {
    result += end_->validate(print);

    count = 0;
    edge = this;
    do
    {
      if( count++ > PST_MAX_LIST_LEN )
      {
        if( print )
        {
          PRINT_ERROR("Bad edge list around Point %p.\n",
                      static_cast<void*>(end_));
        }
        result++;
        break;
      }
      edge = edge->next( end_ );
    } while( edge != this );
  }
  
  if (forward_.face_ && !forward_.next_)
  {
    if( print )
      PRINT_ERROR("Forward CoEdge on Edge %p as Face %p but no next().\n",
                  static_cast<void*>(this), static_cast<void*>(forward_.face_));
    count++;
  }
  
  if (reverse_.face_ && !reverse_.next_)
  {
    if( print )
      PRINT_ERROR("Reverse CoEdge on Edge %p as Face %p but no next().\n",
                  static_cast<void*>(this), static_cast<void*>(reverse_.face_));
    count++;
  }
  
  if (!forward_.face_ && forward_.next_)
  {
    if( print )
      PRINT_ERROR("Forward CoEdge on Edge %p as %p as next CoEdge but no Face.\n",
                  static_cast<void*>(this), static_cast<void*>(forward_.next_));
    count++;
  }
  
  if (!reverse_.face_ && reverse_.next_)
  {
    if( print )
      PRINT_ERROR("Reverse CoEdge on Edge %p as %p as next CoEdge but no Face.\n",
                  static_cast<void*>(this), static_cast<void*>(reverse_.next_));
    count++;
  }
    
  return result;
}

int PST_CoEdge::validate( CubitBoolean print )
{
  int result = 0;
  if( ! edge_ )
  {
    if( print ) PRINT_ERROR("Coedge %p has null edge.\n",
                            static_cast<void*>(this));
    result++;
  }
  else 
  {
    if( ! edge_->other( this ) )
    {
      if( print ) 
        PRINT_ERROR( "Inconsistent coedge->edge link.  CoEdge: %p  Edge: %p\n",
                     static_cast<void*>(edge_), static_cast<void*>(this) );
      result++;
    }
    result += edge_->validate(print);
  }
  return result;
}
      
int PST_Face::validate( CubitBoolean print )
{
  PST_CoEdge* coe = coedge_;
  
  if( !coedge_ )
  {
    if(print) PRINT_ERROR("Face %p has null loop.\n", static_cast<void*>(this));
    return 1;
  }
  
  int count = 0;
  int result = 0;
  do
  {
    if( count++ > PST_MAX_LIST_LEN )
    {
      if( print ) PRINT_ERROR("Face %p loop is infinite.\n",
                              static_cast<void*>(this));
      result++;
      break;
    }
    
    if( coe->face() != this )
    {
      if( print ) 
        PRINT_ERROR("Loop for face %p contains face %p on edge %p.\n",
                    static_cast<void*>(this), static_cast<void*>(coe->face()),
                    static_cast<void*>(coe->edge()) );
      result++;
    }
    
    if( coe->edge() )
    {
      result += coe->edge()->validate( print );
    }
    else
    {
      if( print )
        PRINT_ERROR("Coedge %p in face %p has null edge.\n",
                    static_cast<void*>(coe), static_cast<void*>(this));
      result++;
    }
    
    if( ! coe->next() )
    {
      if( print ) 
        PRINT_ERROR("Null coedge after coedge %p in face %p.\n",
                    static_cast<void*>(coe), static_cast<void*>(this) );
      result++;
      break;
    }

    coe = coe->next();
    
  } while( coe != coedge_ );
  
  return result;
}

int PST_Edge::validate( DLIList<PST_Edge*>& edges, CubitBoolean print )
{
  DLIList<PST_Face*> faces;
  PST_Edge::faces( edges, faces );
  
  int result = 0;
  for( int f = faces.size(); f--; )
    result += faces.get_and_step()->validate(print);
  return result;
}

void PST_Face::print()
{
  if( ! coedge_ )
  {
    PRINT_ERROR("Face %p has null coedge.\n", static_cast<void*>(this));
    return;
  }
  
  PRINT_INFO("Face        CoEdge       Edge         Start Pt     End Pt\n"
             "----------  -----------  -----------  -----------  -----------\n" );
  PST_CoEdge* coe = coedge_;
    PRINT_INFO("%10p %10p  %10p  %10p  %10p\n",
      this, coe, coe->edge(), coe->start_point(), coe->end_point() );
  coe = coe->next();
  while( coe && coe != coedge_ )
  {
    PRINT_INFO("            %10p  %10p  %10p  %10p\n",
      coe, coe->edge(), coe->start_point(), coe->end_point() );
    coe = coe->next();
  }
}

void PST_Edge::print()
{
  PRINT_INFO("Edge %p:  Points:  Start: %p  End: %p   Faces:  Forward: %p  Reverse: %p\n",
             static_cast<void*>(this), static_cast<void*>(start_),
             static_cast<void*>(end_), static_cast<void*>(forward_.face_),
             static_cast<void*>(reverse_.face_) );
}

void PST_Point::print()
{
  PRINT_INFO("Point %p : ( %f, %f, %f )", static_cast<void*>(this), x(), y(),
             z() );
  
  if( edge_ )
  {
    PRINT_INFO("  Edges:");
    int count = 0;
    PST_Edge* edge = edge_;
    do 
    {
      if( count )
      {
        PRINT_INFO("\t%8p", static_cast<void*>(edge));
      }
      else
      {
        PRINT_INFO("\n\t%8p", static_cast<void*>(edge));
      }
      
      count = (count + 1) % 4;
      edge = edge->next(this);
    
    } while( edge && edge != edge_ );
    
  }
  PRINT_INFO("\n");
}

    
PST_Entity::~PST_Entity()
 { }

