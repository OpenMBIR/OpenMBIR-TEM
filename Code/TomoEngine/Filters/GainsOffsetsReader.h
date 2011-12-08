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
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class TomoEngine_EXPORT GainsOffsetsReader : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(GainsOffsetsReader)
    MXA_STATIC_NEW_MACRO(GainsOffsetsReader);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, GainsOffsetsReader);
    MXA_TYPE_MACRO_SUPER(GainsOffsetsReader, TomoFilter)

    virtual ~GainsOffsetsReader();

     virtual void execute();

  protected:
    GainsOffsetsReader();

  private:
    GainsOffsetsReader(const GainsOffsetsReader&); // Copy Constructor Not Implemented
    void operator=(const GainsOffsetsReader&); // Operator '=' Not Implemented
};


#endif /* GAINSOFFSETSREADER_H_ */
