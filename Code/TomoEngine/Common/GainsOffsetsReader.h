/*
 * GainsOffsetsReader.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef GAINSOFFSETSREADER_H_
#define GAINSOFFSETSREADER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class GainsOffsetsReader : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(GainsOffsetsReader)
    MXA_STATIC_NEW_MACRO(GainsOffsetsReader);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, GainsOffsetsReader);
    MXA_TYPE_MACRO_SUPER(GainsOffsetsReader, AbstractFilter)

    virtual ~GainsOffsetsReader();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);



    virtual void execute();

  protected:
    GainsOffsetsReader();

  private:
    GainsOffsetsReader(const GainsOffsetsReader&); // Copy Constructor Not Implemented
    void operator=(const GainsOffsetsReader&); // Operator '=' Not Implemented
};


#endif /* GAINSOFFSETSREADER_H_ */
