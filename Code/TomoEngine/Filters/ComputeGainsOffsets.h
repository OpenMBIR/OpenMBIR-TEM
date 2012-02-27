/*
 * ComputeGainsOffets.h
 *
 *  Created on: Dec 2, 2011
 *      Author: mjackson
 */

#ifndef COMPUTEGAINSOFFETS_H_
#define COMPUTEGAINSOFFETS_H_


#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"
#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"

/*
 *
 */
class TomoEngine_EXPORT ComputeGainsOffsets : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(ComputeGainsOffsets)
    MXA_STATIC_NEW_MACRO(ComputeGainsOffsets);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, ComputeGainsOffsets);
    MXA_TYPE_MACRO_SUPER(ComputeGainsOffsets, TomoFilter)

    virtual ~ComputeGainsOffsets();

    virtual void execute();


  protected:
    ComputeGainsOffsets();

  private:
    ComputeGainsOffsets(const ComputeGainsOffsets&); // Copy Constructor Not Implemented
    void operator=(const ComputeGainsOffsets&); // Operator '=' Not Implemented
};

#endif /* COMPUTEGAINSOFFETS_H_ */
