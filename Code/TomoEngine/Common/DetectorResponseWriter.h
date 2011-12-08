/*
 * DetectorResponseWriterWriter.h
 *
 *  Created on: Dec 8, 2011
 *      Author: mjackson
 */

#ifndef DETECTORRESPONSEWRITER_H_
#define DETECTORRESPONSEWRITER_H_

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"

/*
 *
 */
class TomoEngine_EXPORT DetectorResponseWriter : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(DetectorResponseWriter)
    MXA_STATIC_NEW_MACRO(DetectorResponseWriter);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, DetectorResponseWriter);
    MXA_TYPE_MACRO_SUPER(DetectorResponseWriter, TomoFilter)

    virtual ~DetectorResponseWriter();

    MXA_INSTANCE_PROPERTY(RealVolumeType::Pointer, Response);

    virtual void execute();


  protected:
    DetectorResponseWriter();

  private:
    DetectorResponseWriter(const DetectorResponseWriter&); // Copy Constructor Not Implemented
    void operator=(const DetectorResponseWriter&); // Operator '=' Not Implemented
};

#endif /* DETECTORRESPONSEWRITER_H_ */
