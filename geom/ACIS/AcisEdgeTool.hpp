#ifndef ACIS_EDGE_TOOL_HPP
#define ACIS_EDGE_TOOL_HPP

#include "CubitDefines.h"

class RefEdge;
//class RefVertex;
template <class X> class DLIList;
class AcisModifyEngine;
class Body;
class Curve;
class CubitVector;
class BodySM;
class CurveACIS;

class AcisEdgeTool
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  ~AcisEdgeTool();

  static AcisEdgeTool* instance();
  //- Gives access to the singleton object of this class

  CubitStatus create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr );
  Curve* create_curve_combine( DLIList<Curve*>& curve_list );
  //-  Uses the solid modeller to create a new RefEdge that is a combination 
  //-  of the input chain of edges.  Written for an ACIS method that does this.  
  //-  The input edges are unchanged.

// KGM - remove if really unused.
#if 0
   CubitStatus split_edges_at_vertices( DLIList<RefEdge*>& ref_edge_list, 
                                     DLIList<RefVertex*>& ref_vertex_list,
                                     DLIList<Body*> &new_body_list,
                                     bool keep_old = false );
  CubitStatus split_edges_at_vertices( DLIList<Curve*>& curve_list, 
                                       DLIList<CubitVector*>& point_list,
                                       DLIList<BodySM*> &new_body_list,
                                       bool keep_old = false );
   //-  Splits an edge attached to a body using a vertex location
#endif
  
protected:
   AcisEdgeTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

   static CubitStatus common_setup( const char* name,
                                    DLIList<RefEdge*>& ref_edge_list,
                                    DLIList<Curve*>& curve_list,
                                    DLIList<Body*>& old_bodies );
                                    

private:
  static AcisEdgeTool* instance_;
};

#endif
