/*
 * RawSinogramInitializer.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef RAWSINOGRAMINITIALIZER_H_
#define RAWSINOGRAMINITIALIZER_H_


#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class RawSinogramInitializer : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(RawSinogramInitializer)
    MXA_STATIC_NEW_MACRO(RawSinogramInitializer);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, RawSinogramInitializer);
    MXA_TYPE_MACRO_SUPER(RawSinogramInitializer, AbstractFilter)

    virtual ~RawSinogramInitializer();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);



    virtual void execute();

  protected:
    RawSinogramInitializer();

  private:
    RawSinogramInitializer(const RawSinogramInitializer&); // Copy Constructor Not Implemented
    void operator=(const RawSinogramInitializer&); // Operator '=' Not Implemented
};



#endif /* RAWSINOGRAMINITIALIZER_H_ */
