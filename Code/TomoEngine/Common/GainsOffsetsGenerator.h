/*
 * GainsOffsetsGenerator.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef GAINSOFFSETSGENERATOR_H_
#define GAINSOFFSETSGENERATOR_H_


#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class GainsOffsetsGenerator : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(GainsOffsetsGenerator)
    MXA_STATIC_NEW_MACRO(GainsOffsetsGenerator);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, GainsOffsetsGenerator);
    MXA_TYPE_MACRO_SUPER(GainsOffsetsGenerator, AbstractFilter)

    virtual ~GainsOffsetsGenerator();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);



    virtual void execute();

  protected:
    GainsOffsetsGenerator();

  private:
    GainsOffsetsGenerator(const GainsOffsetsGenerator&); // Copy Constructor Not Implemented
    void operator=(const GainsOffsetsGenerator&); // Operator '=' Not Implemented
};



#endif /* GAINSOFFSETSGENERATOR_H_ */
