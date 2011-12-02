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

#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class ComputeGainsOffsets
{
  public:
    MXA_SHARED_POINTERS(ComputeGainsOffsets)
    MXA_STATIC_NEW_MACRO(ComputeGainsOffsets)
    MXA_TYPE_MACRO(ComputeGainsOffsets);

    virtual ~ComputeGainsOffsets();

    int execute(TomoInputs* inputs, Sinogram* sinogram);


  protected:
    ComputeGainsOffsets();

  private:
    ComputeGainsOffsets(const ComputeGainsOffsets&); // Copy Constructor Not Implemented
    void operator=(const ComputeGainsOffsets&); // Operator '=' Not Implemented
};

#endif /* COMPUTEGAINSOFFETS_H_ */
