/*
 * MRCSinogramInitializer.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef MRCSINOGRAMINITIALIZER_H_
#define MRCSINOGRAMINITIALIZER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class MRCSinogramInitializer : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(MRCSinogramInitializer)
    MXA_STATIC_NEW_MACRO(MRCSinogramInitializer);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, MRCSinogramInitializer);
    MXA_TYPE_MACRO_SUPER(MRCSinogramInitializer, AbstractFilter)

    virtual ~MRCSinogramInitializer();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);



    virtual void execute();

  protected:
    MRCSinogramInitializer();

  private:
    MRCSinogramInitializer(const MRCSinogramInitializer&); // Copy Constructor Not Implemented
    void operator=(const MRCSinogramInitializer&); // Operator '=' Not Implemented
};

#endif /* MRCSINOGRAMINITIALIZER_H_ */
