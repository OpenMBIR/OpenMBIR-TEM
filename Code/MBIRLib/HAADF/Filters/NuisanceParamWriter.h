/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Singanallur Venkatakrishnan (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Pudue
 * Univeristy, BlueQuartz Software nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#ifndef NUISANCEBINWRITER_H_
#define NUISANCEBINWRITER_H_

#include <string>

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/GenericFilters/TomoFilter.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/HAADF/HAADFForwardModel.h"


/**
 * @class NuisanceBinWriter NuisanceBinWriter.h TomoEngine/IO/NuisanceBinWriter.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Dec 9, 2011
 * @version 1.0
 */
class MBIRLib_EXPORT NuisanceParamWriter : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(NuisanceParamWriter)
    MXA_STATIC_NEW_MACRO(NuisanceParamWriter);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, NuisanceParamWriter);
    MXA_TYPE_MACRO_SUPER(NuisanceParamWriter, TomoFilter)

    virtual ~NuisanceParamWriter();

    enum TargetArray {
      Nuisance_I_O,
      Nuisance_mu,
      Nuisance_alpha
    };

    MXA_INSTANCE_PROPERTY(bool, WriteBinary);
    MXA_INSTANCE_STRING_PROPERTY(FileName);
    MXA_INSTANCE_PROPERTY(TargetArray, DataToWrite);
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Data);
    MXA_INSTANCE_PROPERTY(HAADFForwardModel::Pointer, ForwardModel);
    MXA_INSTANCE_PROPERTY(uint16_t, Ntheta);

    void execute();

  protected:
    NuisanceParamWriter();

  private:
    NuisanceParamWriter(const NuisanceParamWriter&); // Copy Constructor Not Implemented
    void operator=(const NuisanceParamWriter&); // Operator '=' Not Implemented
};

#endif /* NUISANCEBINWRITER_H_ */
