/* ============================================================================
 * Copyright (c) 2011, Singanallur Venkatakrishnan <svenkata@purdue.edu>
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
 * Neither the name of Singanallur Venkatakrishnan , Purdue University nor the
 * names of its contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _SCALEOFFSET_CORRECTION_PARSER_H_
#define _SCALEOFFSET_CORRECTION_PARSER_H_

#include <stdio.h>
#include <string>
#include <iostream>


#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/BrightField/BFMultiResolutionReconstruction.h"


/**
 * @class BFReconstructionArgsParser BFReconstructionArgsParser.h ScaleOffsetCorrectionAlgorithm/BFReconstructionArgsParser.h
 * @brief
 * @author
 * @date Nov 16, 2011
 * @version 1.0
 */
class BFReconstructionArgsParser
{
  public:
    BFReconstructionArgsParser();
    virtual ~BFReconstructionArgsParser();

    int parseArguments(int argc, char **argv, BFMultiResolutionReconstruction::Pointer multiRes);

    void printInputs(TomoInputsPtr inputs, std::ostream &out);

  private:
    uint64_t startm;
    uint64_t stopm;

    BFReconstructionArgsParser(const BFReconstructionArgsParser&); // Copy Constructor Not Implemented
    void operator=(const BFReconstructionArgsParser&); // Operator '=' Not Implemented

};

#endif //ComputationInputs
