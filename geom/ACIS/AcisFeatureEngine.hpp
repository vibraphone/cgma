#ifndef _ACISFEATUREENGINE_HPP
#define _ACISFEATUREENGINE_HPP

// *** BEGIN INCLUDES *** //
#include "GeometryFeatureEngine.hpp"
// *** END INCLUDES *** //

// *** BEGIN FORWARD DECLARATIONS *** //
class TopologyBridge;
class TopologyEntity;
class GeometryEntity;
// *** END FORWARD DECLARATIONS *** //

//! Defines an interface for the AcisFeatureEngine
class AcisFeatureEngine : public GeometryFeatureEngine
{
public:
    //! Gives access to the singleton object of this class
    static AcisFeatureEngine* instance();

    //! virtual destructor
    virtual ~AcisFeatureEngine() {}

    //! returns CUBIT_TRUE if the entity belongs to this feature engine
    virtual CubitBoolean is_feature_engine( const TopologyBridge* );

    //! returns CUBIT_TRUE if the entity belongs to this feature engine
    virtual CubitBoolean is_feature_engine( const TopologyEntity* );

    //! Get all of the entities related to 'entity_in' by the same source feature
    virtual CubitStatus get_related_by_feature(
        GeometryEntity *entity_in,
        DLIList<GeometryEntity*> &list_out);

    //! Get the feature types available from this FeatureEngine
    virtual DLIList<GeometryFeatureEngine::FeatureType> 
        available_feature_types(GeometryEntity* entity);

    //! Get the feature type that created 'entity'
    GeometryFeatureEngine::FeatureType feature_type(GeometryEntity* entity);

protected:
    //! protected constructor
    AcisFeatureEngine();

private:
    static AcisFeatureEngine* instance_;
};
#endif
