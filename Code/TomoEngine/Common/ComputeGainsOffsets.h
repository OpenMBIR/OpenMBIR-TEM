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
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"

/*
 *
 */
class ComputeGainsOffsets : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(ComputeGainsOffsets)
    MXA_STATIC_NEW_MACRO(ComputeGainsOffsets)
    MXA_TYPE_MACRO(ComputeGainsOffsets);

    virtual ~ComputeGainsOffsets();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);

    virtual void execute();


  protected:
    ComputeGainsOffsets();

  private:
    ComputeGainsOffsets(const ComputeGainsOffsets&); // Copy Constructor Not Implemented
    void operator=(const ComputeGainsOffsets&); // Operator '=' Not Implemented
};

#endif /* COMPUTEGAINSOFFETS_H_ */
