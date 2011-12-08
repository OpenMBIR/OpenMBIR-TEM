/*
 * DetectorResponse.h
 *
 *  Created on: Dec 8, 2011
 *      Author: mjackson
 */

#ifndef DETECTORRESPONSE_H_
#define DETECTORRESPONSE_H_

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"

/*
 *
 */
class TomoEngine_EXPORT DetectorResponse : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(DetectorResponse)
    MXA_STATIC_NEW_MACRO(DetectorResponse);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, DetectorResponse);
    MXA_TYPE_MACRO_SUPER(DetectorResponse, TomoFilter)

    virtual ~DetectorResponse();

    MXA_INSTANCE_PROPERTY(DATA_TYPE, BeamWidth);
    MXA_INSTANCE_PROPERTY(DATA_TYPE, OffsetR);
    MXA_INSTANCE_PROPERTY(DATA_TYPE, OffsetT);
    MXA_INSTANCE_PROPERTY(RealImageType::Pointer, VoxelProfile);
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, BeamProfile);
    MXA_INSTANCE_PROPERTY(RealVolumeType::Pointer, Response);

    virtual void execute();


  protected:
    DetectorResponse();

  private:
    DetectorResponse(const DetectorResponse&); // Copy Constructor Not Implemented
    void operator=(const DetectorResponse&); // Operator '=' Not Implemented
};


#endif /* DETECTORRESPONSE_H_ */
