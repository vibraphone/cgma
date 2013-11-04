#ifndef INTERMEDIATE_GEOM_ENGINE_HPP
#define INTERMEDIATE_GEOM_ENGINE_HPP

template <class X> class DLIList;
class TopologyBridge;
class CubitTransformMatrix;
class Body;
class Surface;
class Curve;
class TBPoint;
class BodySM;
class TBOwner;

class IntermediateGeomEngine
{

public:
  virtual bool is_composite(TBOwner *bridge_owner) = 0;
  virtual bool is_composite(TopologyBridge *bridge ) = 0;
  virtual bool is_partition(TBOwner *bridge_owner) = 0;

  virtual int level() const = 0;
  
  virtual void remove_imprint_attributes_after_modify
                                ( DLIList<BodySM*> &old_sms,
                                  DLIList<BodySM*> &new_sms )=0;
  virtual void push_imprint_attributes_before_modify
                                ( DLIList<BodySM*> &body_sms ) = 0;
  virtual void push_named_attributes_to_curves_and_points
                                ( DLIList<TopologyBridge*> &tb_list, const char *name_in ) = 0;
  virtual CubitStatus export_geometry( DLIList<TopologyBridge*>& geometry_list ) = 0;
  
  virtual CubitStatus import_geometry( DLIList<TopologyBridge*>& geometry_list ) = 0;
  
  virtual void clean_out_deactivated_geometry() = 0;

  virtual void remove_attributes( DLIList<TopologyBridge*> &bridge_list ) = 0;
  virtual void attribute_after_imprinting(DLIList<TopologyBridge*> &tb_list,
                                          DLIList<Body*> &old_bodies)=0;
  
  virtual void remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges) = 0;
  virtual void remove_modified(DLIList<Surface*> &all_surfs,
    DLIList<Curve*> &all_curves, DLIList<TBPoint*> &all_pts) = 0;
  virtual CubitStatus notify_transform( TopologyBridge* entity,
                                        const CubitTransformMatrix& xform ) = 0;

  virtual void get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                               DLIList<TopologyBridge*> &tbs ) = 0;

};

#endif

