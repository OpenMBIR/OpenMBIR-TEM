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
#include "TomoEngine/Common/InitialReconstructionInitializer.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class TomoEngine_EXPORT InitialReconstructionBinReader : public InitialReconstructionInitializer
{
  public:
    MXA_SHARED_POINTERS(InitialReconstructionBinReader)
    MXA_STATIC_NEW_MACRO(InitialReconstructionBinReader);
    MXA_STATIC_NEW_SUPERCLASS(InitialReconstructionInitializer, InitialReconstructionBinReader);
    MXA_TYPE_MACRO_SUPER(InitialReconstructionBinReader, InitialReconstructionInitializer)

    virtual ~InitialReconstructionBinReader();



    virtual void execute();

    virtual void initializeData();

  protected:
    InitialReconstructionBinReader();

  private:
    InitialReconstructionBinReader(const InitialReconstructionBinReader&); // Copy Constructor Not Implemented
    void operator=(const InitialReconstructionBinReader&); // Operator '=' Not Implemented
};


#endif /* INITIALRECONSTRUCTIONBINREADER_H_ */
