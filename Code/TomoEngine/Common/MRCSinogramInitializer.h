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
#include "TomoEngine/Common/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class TomoEngine_EXPORT MRCSinogramInitializer : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(MRCSinogramInitializer)
    MXA_STATIC_NEW_MACRO(MRCSinogramInitializer);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, MRCSinogramInitializer);
    MXA_TYPE_MACRO_SUPER(MRCSinogramInitializer, TomoFilter)

    virtual ~MRCSinogramInitializer();

    virtual void execute();

  protected:
    MRCSinogramInitializer();

  private:
    MRCSinogramInitializer(const MRCSinogramInitializer&); // Copy Constructor Not Implemented
    void operator=(const MRCSinogramInitializer&); // Operator '=' Not Implemented
};

#endif /* MRCSINOGRAMINITIALIZER_H_ */
