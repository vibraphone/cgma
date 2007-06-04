#ifndef ACIS_DRAW_TOOL_HPP
#define ACIS_DRAW_TOOL_HPP

#include "CubitDefines.h"

class ENTITY;
class FACE;
class EDGE;
class VERTEX;
class Surface;
class Curve;

template <class X> class DLIList;

class AcisDrawTool
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  ~AcisDrawTool();

  static AcisDrawTool* instance();
  //- Gives access to the singleton object of this class

  //do not remove this function....it is handy for debugging.
  CubitStatus draw_ENTITY( ENTITY *ENTITY_ptr, int color,
                           CubitBoolean tessellate = CUBIT_FALSE,
                           CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given ENTITY.  Handles entities of type BODY, LUMP, SHELL,
  //- FACE, LOOP, COEDGE, EDGE or VERTEX.  The tessellate option only applies
  //- to FACEs - if true, the FACE is tessellated and each triangle is drawn;
  //- otherwise, only the outer curves are drawn.

  CubitStatus draw_surface( Surface *surface, int color, CubitBoolean tessellate );
  CubitStatus draw_curve( Curve *curve, int color );

  CubitStatus draw_FACE( FACE *FACE_ptr, int color, 
                         CubitBoolean tessellate = CUBIT_FALSE,
                         CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given FACE.  If tessellate is true, the FACE is tessellated and
  //- each triangle is drawn; otherwise, only the outer curves are drawn.

  CubitStatus draw_EDGE( EDGE *EDGE_ptr, int color, 
                         CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given EDGE

  CubitStatus draw_VERTEX( VERTEX *VERTEX_ptr, int color, 
                           CubitBoolean flush = CUBIT_FALSE );
  //- Draw the given VERTEX

protected:
   AcisDrawTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

private:
  
  static AcisDrawTool* instance_;
};

#endif
