/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Singanallur Venkatakrishnan (Purdue University)
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
#ifndef NUISANCEPARAMREADER_H_
#define NUISANCEPARAMREADER_H_

#include <string>

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"

#include "ReconstructionCoreLib/ReconstructionCoreLib.h"
#include "ReconstructionCoreLib/Common/TomoArray.hpp"
#include "HaadfMbirLib/HaadfMbirStructures.h"
#include "ReconstructionCoreLib/GenericFilters/TomoFilter.h"


/*
 *
 */
class ReconstructionCoreLib_EXPORT NuisanceParamReader : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(NuisanceParamReader)
    MXA_STATIC_NEW_MACRO(NuisanceParamReader);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, NuisanceParamReader);
    MXA_TYPE_MACRO_SUPER(NuisanceParamReader, TomoFilter)

    virtual ~NuisanceParamReader();

//    enum TargetArray {
//      Nuisance_I_O,
//      Nuisance_mu,
//      Nuisance_alpha
//    };

  //  MXA_INSTANCE_PROPERTY(bool, WriteBinary);
    MXA_INSTANCE_STRING_PROPERTY(FileName);
  //  MXA_INSTANCE_PROPERTY(TargetArray, DataToRead);
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Data);

    void execute();

  protected:
    NuisanceParamReader();

  private:
    NuisanceParamReader(const NuisanceParamReader&); // Copy Constructor Not Implemented
    void operator=(const NuisanceParamReader&); // Operator '=' Not Implemented

};

#endif /* NUISANCEPARAMREADER_H_ */
