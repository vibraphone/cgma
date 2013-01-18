//-------------------------------------------------------------------------
// Filename      : PST_Data.hpp
//
// Purpose       : Classes for representing geometry facets.  Used by 
//                 PartSurfTess to find the graphis faceting for 
//                 PartitionSurfaces.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/17/01
//-------------------------------------------------------------------------
#ifndef PST_DATA_HPP
#define PST_DATA_HPP
#include <vector>

#include "CubitVector.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitPlane.hpp"

#define PST_OTHER(A,B,C) (((A)==(B))?(C):(((A)==(C))?(B):0))

class PST_CoEdge;
class PST_Edge;
class PST_Point;
class PST_Face;
class GMem;

class PST_Entity;

class PST_EntityOwner {
  public:
  virtual void 
  notify_split( PST_Entity* old_entity, PST_Entity* new_entity ) = 0;
};

class PST_Entity
{
  public:
  
  inline PST_EntityOwner* owner() const;
  
  inline void owner( PST_EntityOwner* );
  
  protected:
  
  inline PST_Entity();
  
  virtual ~PST_Entity() = 0;
  
  private:
  
  PST_EntityOwner* myOwner;
};


//-------------------------------------------------------------------------
// Purpose       : Class representing a vertex in facet data structure
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/01
//-------------------------------------------------------------------------
class PST_Point : public CubitVector, public PST_Entity
{
  friend class PST_Edge;

  public:
  
    PST_Point( const CubitVector& position )
      : CubitVector( position ),
        mark(0), private_mark_(false), edge_(0)
      { }

    ~PST_Point() ;
      //- Destroys this point, and all edges and faces that
      //- the point exists in.
    
    PST_Edge* edge() 
      { return edge_; }
    
    inline PST_Edge* next( PST_Edge* prev );
    
    PST_Edge* common( PST_Point* pt );
    
    int faces( DLIList<PST_Face*>* list = 0 );
    
    //PST_Face* boundary_face();
    
    void debug_draw( int color = 1, bool flush = true );
    static void debug_draw_points( DLIList<PST_Point*>& list, 
                             int color, bool flush = true );
  
    int validate( CubitBoolean print = CUBIT_FALSE );
    
    void print();

    int mark;
  
    int sequence;
    
  private:
    cBit private_mark_ : 1;
  
    PST_Edge* edge_;
};

//-------------------------------------------------------------------------
// Purpose       : Association between a PST_Edge and a PST_Face
//
// Special Notes : Allocated inline in PST_Edge -- Do NOT memory manage!
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/01
//-------------------------------------------------------------------------
class PST_CoEdge
{
  friend class PST_Edge;
  
  public:
  
    PST_CoEdge* next()
      { return next_; }
      
    PST_Edge* edge()
      { return edge_; }
    
    PST_Face* face()
      { return face_; }
    
    inline CubitSense sense();
    
    inline PST_Point* start_point();
    inline PST_Point* end_point();
    inline PST_Point* common( PST_CoEdge* coedge );
    inline PST_Point* other( PST_Point* point );
    
    inline const CubitVector& start_coord();
    inline const CubitVector& end_coord();
    inline CubitVector direction();
    
    inline PST_CoEdge* other();
    
    PST_CoEdge* previous();
      
    int validate( CubitBoolean print = CUBIT_FALSE );
      
  private:
  
    PST_CoEdge( PST_Edge* owner )
      : next_(0), edge_(owner), face_(0)
      { }
    
    ~PST_CoEdge() 
      { assert( !face_ ); } 
    
    PST_CoEdge* next_;
    PST_Edge*   edge_;
    PST_Face*   face_;
};
    

//-------------------------------------------------------------------------
// Purpose       : Edge in facet data structure
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/01
//-------------------------------------------------------------------------
class PST_Edge : public PST_Entity
{
  public:
    
    static void make_gmem( GMem& gmem, DLIList<PST_Face*>& facets );
    static void make_gmem( GMem& gmem, DLIList<PST_Edge*>& edges );
    static void make_facets( GMem& gmem, double tolerance,
                             DLIList<PST_Edge*>& edges );
    static void make_facets(
     			     const std::vector<double>& coordinates, 
               const std::vector<int>& connections,
    			     double tolerance, 
               DLIList<PST_Edge*>& edge_list );

    PST_Edge( PST_Point* start_pt, PST_Point* end_pt )
      : mark(0), private_mark_(false),
        start_(0), end_(0),
        forward_( this ), reverse_( this )
      { 
        assert( start_pt != end_pt );
        set_start_point( start_pt );
        set_end_point( end_pt );
      }
    
    
    static PST_Face* create_face( PST_Point* p1,
                                  PST_Point* p2,
                                  PST_Point* p3 );

    static PST_Edge* split_face( PST_Point* start,
                                 PST_Point* end,
                                 PST_Face* hint = 0 );
      //- This method splits a face by creating an edge
      //- between two existing points on that face.   
      //- If the face to split is not specified, both
      //- points must have exactly one face in common.
      //- If the points do not have a common face, have
      //- more than one common face and no face is 
      //- specified, do not exist in the passed face, 
      //- or already have a common edge, NULL is returned.
                                 
    static PST_Edge* split_face( PST_CoEdge* after_coedge_1,
                                 PST_CoEdge* after_coedge_2 );
      //- Split a face with an edge from the end of the first
      //- coedge to the end of the second coedge.  NULL is
      //- returned if the coedges do not have a common face,
      //- or the desired edge already exists.
                                 
    static PST_Edge* insert_in_face( PST_Point* end,
                                     PST_CoEdge* after_this );
      //- Create a non-manifold edge in a face.  The edge is
      //- inserted starting at the end of the passed coedge,
      //- and terminating at the passed point.  If the 
      //- passed point already exists in a face, null is
      //- returned.
                                     
    PST_Edge* split( PST_Point* point );
      //- Split an existing edge.  Creates a new edge in the 
      //- same face(s) as the existing edge.  Returns null
      //- if the passed already exists in any face.

    ~PST_Edge();
      //- Destroy this edge, splitting or merging faces as
      //- necessary.  Does NOT destroy unused points.
    
    PST_Point* start_point()
      { return start_; }
    
    PST_Point* end_point()
      { return end_; }
    
    const CubitVector& start_coord()
      { return *start_; }
      
    const CubitVector& end_coord()
      { return *end_; }
    
    CubitVector direction()
      { return end_coord() - start_coord(); }
      
    
    CubitVector position( double t )
      { return start_coord() + t * direction(); }
  
    double closest_on_line( const CubitVector& pos );
    double closest_on_edge( const CubitVector& pos );
      // Parameter of the closest position on 
      // a) the line containing this edge or 
      // b) this bounded edge.
      
    double closest_on_line( const CubitVector& base, 
                            const CubitVector& direction );
      // The parameter of the closest position on
      // the line containing this edge to the passed line.
      // Use position(..) function above to get coordinates.
      
    PST_Point* other( PST_Point* pt )
      { return PST_OTHER(pt,start_,end_); }
    
    inline PST_Point* common_point( PST_Edge* edge );
    
    PST_CoEdge* forward() 
      { return &forward_; }
    
    PST_CoEdge* reverse()
      { return &reverse_; }
      
    PST_CoEdge* other( PST_CoEdge* coedge )
      { return PST_OTHER(coedge,&forward_,&reverse_); }
    
    inline PST_CoEdge* coedge( PST_Face* face );
    
    
    PST_Face* other( PST_Face* face )
      { return PST_OTHER(face,forward_.face(),reverse_.face()); }
   
    
    inline PST_Face* common_face( PST_Edge* edge );
    
    inline CubitSense sense( PST_Face* face );
    
    CubitSense sense( PST_CoEdge* coedge )
      { return coedge == forward() ? CUBIT_FORWARD  :
               coedge == reverse() ? CUBIT_REVERSED : CUBIT_UNKNOWN; }

    PST_Edge* other( PST_Point* point, PST_Face* face );
      // Find the other edge in the passed face that
      // shares the passed point.  Returns NULL if multiple
      // matches in a non-manifold face.

    PST_Edge* next( PST_Point* pt )
      { return pt == start_ ? start_next_ : 
               pt == end_   ? end_next_ : 0; }
               
    //inline bool boundary();
    
    static void faces( DLIList<PST_Edge*>& edges, DLIList<PST_Face*>& faces );
    static void edges( DLIList<PST_Face*>& faces, DLIList<PST_Edge*>& edges );
    static void edges( DLIList<PST_Point*>& pts,  DLIList<PST_Edge*>& edges );
    static void points( DLIList<PST_Edge*>& edges, DLIList<PST_Point*>& pts );
    
    void debug_draw( int color = 1, bool flush = true );
    static void debug_draw_points( DLIList<PST_Edge*>& list, 
                             int color, int bcolor = 0, bool flush = true );
    static void debug_draw_edges ( DLIList<PST_Edge*>& list, 
                             int color, int bcolor = 0, bool flush = true );
    static void debug_draw_faces ( DLIList<PST_Edge*>& list, 
                             int color, bool flush = true );
    
    static int validate( DLIList<PST_Edge*>& edges, CubitBoolean print = CUBIT_FALSE );
    int validate( CubitBoolean print = CUBIT_FALSE );
    
    void print();

    int mark;

  private:
    cBit private_mark_ : 1;
    
    void remove_point( PST_Point*& ptr, PST_Edge*& next );
    void remove_start_point() { remove_point( start_, start_next_ ); }
    void remove_end_point()   { remove_point( end_, end_next_); }
    
    void set_point( PST_Point* pt, PST_Point*& ptr, PST_Edge*& next );
    void set_start_point( PST_Point* pt )
      { set_point( pt, start_, start_next_ ); }
    void set_end_point( PST_Point* pt )
      { set_point( pt, end_, end_next_ ); }
      
    PST_Point* start_;
    PST_Point* end_;
    PST_Edge* start_next_;
    PST_Edge* end_next_;
    
    PST_CoEdge forward_;
    PST_CoEdge reverse_;
};

//-------------------------------------------------------------------------
// Purpose       : Face in facet data structure
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/01
//-------------------------------------------------------------------------
class PST_Face : public PST_Entity
{
  friend class PST_Edge;
  friend class PST_Point;

  public:
      
    const CubitPlane& plane()
      {
        if( update_plane_ )
          calculate_plane();
        return plane_;
      }
  
    CubitVector normal()
      { return plane().normal(); }
      
      
    double bounding_length();
    
    PST_CoEdge* first_coedge()
      { return coedge_; }
      
    void append_points( DLIList<PST_Point*>& result_list );
    
    PST_Point* opposite( PST_Edge* edge );
    PST_Edge* opposite( PST_Point* point );
      // for triangles
      
    bool two_edges( PST_Point* point, PST_Edge*& edge1, PST_Edge*& edge2 );
  
    bool modified()
      { return modified_; }
    
    void modified( bool b )
      { 
        modified_ = b;
        if( b ) update_plane_ = 1;
      }
      
    void debug_draw( int color = 1, bool flush = true );
    static void debug_draw_faces( DLIList<PST_Face*>& list, 
                            int color, bool flush = true );
  
    int validate( CubitBoolean print = CUBIT_FALSE );
    
    void print();
    
    int mark;
    
    int sequence;
    
    int parent;
    
  private:
  
    PST_Face( PST_CoEdge* coedge )
      : mark(0), sequence(0), parent(0),
        coedge_(coedge), 
        //boundary_(0),
        modified_(0),
        update_plane_(1),
        private_mark_(false)
        
      { }
      
    PST_Face( PST_CoEdge* coedge, PST_Face* split_from )
      : mark(0), sequence(0), parent(0),
        coedge_(coedge),
        //boundary_(split_from->boundary()),
        modified_(0),
        update_plane_(1),
        private_mark_(false)
      { split_from->modified(true); }
  
    ~PST_Face()
      { assert( !coedge_ ); }
  
    bool calculate_plane();
    
    PST_CoEdge* coedge_;
    CubitPlane  plane_;
    
    //cBit boundary_     : 1;
    cBit modified_     : 1;
    cBit update_plane_ : 1;
    cBit private_mark_ : 1;
};

/*************************************************************************
 *                               PST_Point
 ************************************************************************/

inline PST_Edge* PST_Point::next( PST_Edge* prev )
{ 
  return prev->next(this);
}

/*************************************************************************
 *                               PST_CoEdge
 ************************************************************************/

inline CubitSense PST_CoEdge::sense()
{
  return edge_->sense(this);
}

inline PST_Point* PST_CoEdge::start_point()
{
  return sense() == CUBIT_FORWARD ? edge_->start_point() : edge_->end_point();
}

inline PST_Point* PST_CoEdge::end_point()
{
  return sense() == CUBIT_FORWARD ? edge_->end_point() : edge_->start_point();
}

inline PST_Point* PST_CoEdge::other( PST_Point* point )
{
  return edge_->other( point );
}

inline const CubitVector& PST_CoEdge::start_coord()
{
  return sense() == CUBIT_FORWARD ? edge_->start_coord() : edge_->end_coord();
}

inline const CubitVector& PST_CoEdge::end_coord()
{
  return sense() == CUBIT_FORWARD ? edge_->end_coord() : edge_->start_coord();
}

inline CubitVector PST_CoEdge::direction()
{
  return sense() == CUBIT_FORWARD                  ? 
         edge_->end_coord() - edge_->start_coord() :
         edge_->start_coord() - edge_->end_coord() ;
}

inline PST_CoEdge* PST_CoEdge::other()
{
  return edge_->other(this);
}



/*************************************************************************
 *                               PST_Edge
 ************************************************************************/

inline PST_Point* PST_Edge::common_point( PST_Edge* edge )
{
  PST_Point* result = 0;
  
  if( (edge->start_ == start_) || 
      (edge->end_   == start_)  )
    result = start_;
  
  if( (edge->start_ == end_) || 
      (edge->end_   == end_)  )
    result = result ? 0 : end_;
  
  return result;
}

inline PST_CoEdge* PST_Edge::coedge( PST_Face* face )
{
  PST_CoEdge* result = 0;
  
  if( forward_.face_ == face )
    result = &forward_;
  
  if( reverse_.face_ == face )
    result = result ? 0 : &reverse_;
  
  return result;
}

inline CubitSense PST_Edge::sense( PST_Face* face )
{
  CubitSense result = CUBIT_UNKNOWN;
  
  if( forward_.face_ == face )
    result = CUBIT_FORWARD;
  
  if( reverse_.face_ == face )
    result = (result != CUBIT_UNKNOWN) ? CUBIT_UNKNOWN : CUBIT_REVERSED;
  
  return result;
}

inline PST_Entity::PST_Entity() : myOwner(0) { ; }

inline PST_EntityOwner* PST_Entity::owner() const
  { return myOwner; }

inline void PST_Entity::owner( PST_EntityOwner* ent )
  { myOwner = ent; }


/*************************************************************************
 *                               PST_Face
 ************************************************************************/
 
 // find the minimum edge length and use this for the bounding length
inline double PST_Face::bounding_length()
{
  // The coedges appear to be stored in a circular list so 
  // get the length of the first edge and then iterate
  PST_CoEdge* coe = this->first_coedge();
  CubitVector pos1 = coe->start_coord();
  CubitVector pos2 = coe->end_coord();
  CubitVector dist = pos1 - pos2;
  double length = dist.length_squared();

  // start the chain at the next coedge
  coe = coe->next();

  // walk through the coedges till we get back to the first 
  // and find the minimum length edge
  double l;
  for ( ; coe != this->first_coedge(); coe = coe->next() )
  {
    pos1 = coe->start_coord();
    pos2 = coe->end_coord();
    dist = pos1 - pos2;
    l = dist.length_squared();
    if ( l < length)
      length = l;
  }

  return sqrt(length);
}

#endif

