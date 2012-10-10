//---------------------------------------------------------------------------
//- Filename:       CollapseCurveTool
//- Purpose:  To collapse small curves for preparing for mesh
//-
//- Creator:       Brett Clark
//- Creation date: 01/10/2006
//---------------------------------------------------------------------------- 

#ifndef COLLAPSECURVETOOL_HPP
#define COLLAPSECURVETOOL_HPP

class RefEdge;
class RefVertex;
template<class X> class DLIList;
   
class CollapseCurveTool
{
  public:
    static CollapseCurveTool *instance();
    static void delete_instance() { if(instance_) delete instance_; };

    CubitStatus collapse_curve(DLIList <RefEdge*> ref_edge_list, 
                DLIList<RefVertex*> ref_vertex_list,
                int ignore_surfaces,
                double *small_curve_size = NULL); 

    CubitStatus collapse_curve_only_virtual(DLIList <RefEdge*> ref_edge_list, 
                                        DLIList<RefVertex*> ref_vertex_list,
                                        int ignore_surfaces,
                                        double *small_curve_size = NULL);
  private:
    CollapseCurveTool();
    // Constructor
                                                                                
    ~CollapseCurveTool();
    // Destructor
                                                                                
    static CollapseCurveTool *instance_;
    // the static instance pointer

    CubitStatus position_from_length(RefEdge *edge,
                                     RefVertex *root_vertex,
                                     double  arc_length,
                                     CubitVector& v_new);


};

#endif

