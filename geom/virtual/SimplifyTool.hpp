//-------------------------------------------------------------------------
// Filename      : SimplifyTool.hpp
//
// Purpose       :
//
// Simplify {Volume|Surface} <Range> [Angle <Value>]
//     [Respect {Surface <Range> | Curve <Range> | Fillet |Chamfer | imprint}]
//     [Preview]
//
// Special Notes :
//
// Creator       : Sam Showman
//
// Creation Date : 11/07/2005
//-------------------------------------------------------------------------
#ifndef SIMPLIFYTOOL_HPP
#define SIMPLIFYTOOL_HPP

#include "CubitDefines.h"
#include "CubitVector.hpp"

class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
template <class X> class DLIList;

/// Tool to automatically simplify geometry
/**
*/
class SimplifyTool
{
public:

    SimplifyTool();
    ~SimplifyTool();

    /// Simplify a list of volumes
    /***/
    CubitStatus simplify_volumes(
        DLIList<RefVolume*> ref_volume_list,
        double surf_angle_in,
        DLIList<RefFace*> respect_face_list,
        DLIList<RefEdge*> respect_edge_list,
        CubitBoolean respect_rounds,
        CubitBoolean respect_imprints,
        CubitBoolean local_normals,
        CubitBoolean preview);

    CubitStatus simplify_surfaces(
        DLIList<RefFace*> ref_face_list, 
        double angle_in,
        DLIList<RefFace*> respect_face_list,
        DLIList<RefEdge*> respect_edge_list,
        CubitBoolean respect_rounds,
        CubitBoolean respect_imprints,
        CubitBoolean local_normals,
        CubitBoolean preview);

	CubitStatus simplify_curves(
		DLIList<RefEdge*> ref_edge_list, 
		double angle_in,
		DLIList<RefEdge*> respect_edge_list,
		DLIList<RefVertex*> respect_vertex_list,
		CubitBoolean respect_imprints,
                    CubitBoolean local_normals,
		CubitBoolean preview);

	CubitStatus simplify_curves_in_volume(
		DLIList<RefEdge*> ref_edge_list, 
		double angle_in,
		DLIList<RefEdge*> respect_edge_list,
		DLIList<RefVertex*> respect_vertex_list,
		CubitBoolean respect_imprints,
                    CubitBoolean local_normals,
		CubitBoolean preview);

	CubitBoolean composite_curves(
		RefEdge* seed_ref_edge,
		RefEdge* ref_edge_ptr,
		double angle_in);

private:

    CubitStatus simplify_volume(
        RefVolume* ref_volume_list, 
        double surf_angle_in,
        DLIList<RefFace*> respect_face_list,
        DLIList<RefEdge*> respect_edge_list,
        CubitBoolean respect_rounds,
        CubitBoolean respect_imprints,
        CubitBoolean local_normals,
        CubitBoolean preview);

    CubitStatus simplify_surfaces_in_volume(
        DLIList<RefFace*> ref_face_list, 
        double angle_in,
        DLIList<RefFace*> respect_face_list,
        DLIList<RefEdge*> respect_edge_list,
        CubitBoolean respect_rounds,
        CubitBoolean respect_imprints,
        CubitBoolean local_normals,
        CubitBoolean preview);

    CubitBoolean composite_surfaces(
        RefFace* seed_ref_face,
        RefFace* ref_face_ptr,
        double angle_in);

    void process_rounds(
        RefVolume* ref_volume,
        double min_radius = 0.0, 
        double max_radius = 1.0);

    CubitStatus
        weighted_average_normal(RefFace* ref_face,
        CubitVector &normal, 
        double &weight );

  
  CubitBoolean composite_surfaces_test_at_curves(RefFace* seed_ref_face,
                                                RefFace* ref_face_ptr,
                                                double angle_in);
  


  CubitBoolean maximum_angle_deviation(RefFace* seed_ref_face,
                                       RefFace* ref_face_ptr,
                                       double angle_in,
                                       double &angle_out);
  

};

#endif

