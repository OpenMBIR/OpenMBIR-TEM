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

#ifndef INITIALRECONSTRUCTIONINITIALIZER_H_
#define INITIALRECONSTRUCTIONINITIALIZER_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "ReconstructionCoreLib/ReconstructionCoreLib.h"
#include "ReconstructionCoreLib/GenericFilters/TomoFilter.h"
#include "HaadfMbirLib/HaadfMbirStructures.h"


/**
 * @class InitialReconstructionInitializer InitialReconstructionInitializer.h SOC/InitialReconstructionInitializer.h
 * @brief
 * @author
 * @date
 * @version 1.0
 */
class ReconstructionCoreLib_EXPORT InitialReconstructionInitializer : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(InitialReconstructionInitializer)
    MXA_STATIC_NEW_MACRO(InitialReconstructionInitializer);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, InitialReconstructionInitializer);
    MXA_TYPE_MACRO_SUPER(InitialReconstructionInitializer, TomoFilter)

    virtual ~InitialReconstructionInitializer();

    Real_t absMaxArray(std::vector<Real_t> &Array);

    virtual void execute();

    virtual void initializeData();

  protected:
    InitialReconstructionInitializer();

  private:
    InitialReconstructionInitializer(const InitialReconstructionInitializer&); // Copy Constructor Not Implemented
    void operator=(const InitialReconstructionInitializer&); // Operator '=' Not Implemented
};


#endif /* INITIALRECONSTRUCTIONINITIALIZER_H_ */
