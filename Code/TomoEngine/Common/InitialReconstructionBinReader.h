/*
 * InitialReconstructionBinReader.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef INITIALRECONSTRUCTIONBINREADER_H_
#define INITIALRECONSTRUCTIONBINREADER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class InitialReconstructionBinReader : public AbstractFilter
{
  public:
    MXA_SHARED_POINTERS(InitialReconstructionBinReader)
    MXA_STATIC_NEW_MACRO(InitialReconstructionBinReader);
    MXA_STATIC_NEW_SUPERCLASS(AbstractFilter, InitialReconstructionBinReader);
    MXA_TYPE_MACRO_SUPER(InitialReconstructionBinReader, AbstractFilter)

    virtual ~InitialReconstructionBinReader();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);


    DATA_TYPE absMaxArray(std::vector<DATA_TYPE> &Array);


    virtual void execute();

  protected:
    InitialReconstructionBinReader();

  private:
    InitialReconstructionBinReader(const InitialReconstructionBinReader&); // Copy Constructor Not Implemented
    void operator=(const InitialReconstructionBinReader&); // Operator '=' Not Implemented
};


#endif /* INITIALRECONSTRUCTIONBINREADER_H_ */
