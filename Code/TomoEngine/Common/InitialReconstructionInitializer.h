/*
 * InitialReconstructionInitializer.h
 *
 *  Created on: Dec 6, 2011
 *      Author: mjackson
 */

#ifndef INITIALRECONSTRUCTIONINITIALIZER_H_
#define INITIALRECONSTRUCTIONINITIALIZER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class InitialReconstructionInitializer : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(InitialReconstructionInitializer)
    MXA_STATIC_NEW_MACRO(InitialReconstructionInitializer);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, InitialReconstructionInitializer);
    MXA_TYPE_MACRO_SUPER(InitialReconstructionInitializer, TomoFilter)

    virtual ~InitialReconstructionInitializer();

    DATA_TYPE absMaxArray(std::vector<DATA_TYPE> &Array);

    virtual void execute();

    virtual void initializeData();

  protected:
    InitialReconstructionInitializer();

  private:
    InitialReconstructionInitializer(const InitialReconstructionInitializer&); // Copy Constructor Not Implemented
    void operator=(const InitialReconstructionInitializer&); // Operator '=' Not Implemented
};


#endif /* INITIALRECONSTRUCTIONINITIALIZER_H_ */
