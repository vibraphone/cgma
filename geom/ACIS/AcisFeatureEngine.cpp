// ********** BEGIN STANDARD INCLUDES **********
// ********** END STANDARD INCLUDES   **********

// ********** BEGIN ACIS Includes  **********
#include "AcisFeatureEngine.hpp"
// ********** END ACIS Includes    **********

// ********** BEGIN CUBIT INCLUDES    **********
#include "TopologyBridge.hpp"
#include "TopologyEntity.hpp"
#include "GeometryFeatureTool.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "AcisBridge.hpp"
#include "CurveACIS.hpp"
// ********** END CUBIT INCLUDES      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

AcisFeatureEngine* AcisFeatureEngine::instance_ = 0;

AcisFeatureEngine* AcisFeatureEngine::instance()
{
    if( instance_ == 0 )
        new AcisFeatureEngine();

    return instance_;
}

AcisFeatureEngine::AcisFeatureEngine()
{
    instance_ = this;
    GeometryFeatureTool::instance()->add_gfe(this);
}

CubitBoolean AcisFeatureEngine::is_feature_engine( const TopologyBridge* tb_ptr )
{
    if(dynamic_cast<AcisBridge*>(const_cast<TopologyBridge*>(tb_ptr))) 
        return CUBIT_TRUE;
    else 
        return CUBIT_FALSE;
}

CubitBoolean AcisFeatureEngine::is_feature_engine( const TopologyEntity* te_ptr )
{
    TopologyBridge *tb_ptr = te_ptr->bridge_manager()->topology_bridge();
    return is_feature_engine(tb_ptr);
}

CubitStatus AcisFeatureEngine::get_related_by_feature(GeometryEntity *entity_in,
                                                      DLIList<GeometryEntity*> &list_out)
{
    return CUBIT_FAILURE;
}

DLIList<GeometryFeatureEngine::FeatureType> 
AcisFeatureEngine::available_feature_types(GeometryEntity* entity)
{
    DLIList<FeatureType> feature_info_list;
    feature_info_list.append(FEATURE_IMPRINT);
    return feature_info_list;
}

GeometryFeatureEngine::FeatureType
AcisFeatureEngine::feature_type(GeometryEntity* entity)
{
    // Because I'm using TDSourceFeature I'm not looking for the feature here
    return FEATURE_UNDEFINED;
}
