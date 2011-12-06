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
#include "TomoEngine/Common/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class RawSinogramInitializer : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(RawSinogramInitializer)
    MXA_STATIC_NEW_MACRO(RawSinogramInitializer);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, RawSinogramInitializer);
    MXA_TYPE_MACRO_SUPER(RawSinogramInitializer, TomoFilter)

    virtual ~RawSinogramInitializer();

    virtual void execute();

  protected:
    RawSinogramInitializer();

  private:
    RawSinogramInitializer(const RawSinogramInitializer&); // Copy Constructor Not Implemented
    void operator=(const RawSinogramInitializer&); // Operator '=' Not Implemented
};



#endif /* RAWSINOGRAMINITIALIZER_H_ */
