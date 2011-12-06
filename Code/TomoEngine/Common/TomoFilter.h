/*
 * TomoFilter.h
 *
 *  Created on: Dec 6, 2011
 *      Author: mjackson
 */

#ifndef TOMOFILTER_H_
#define TOMOFILTER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class TomoFilter : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(TomoFilter)
    MXA_STATIC_NEW_MACRO(TomoFilter);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, TomoFilter);
    MXA_TYPE_MACRO_SUPER(TomoFilter, AbstractFilter)

    virtual ~TomoFilter();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);

    virtual void execute();

  protected:
    TomoFilter();

  private:
    TomoFilter(const TomoFilter&); // Copy Constructor Not Implemented
    void operator=(const TomoFilter&); // Operator '=' Not Implemented
};


#endif /* TOMOFILTER_H_ */
