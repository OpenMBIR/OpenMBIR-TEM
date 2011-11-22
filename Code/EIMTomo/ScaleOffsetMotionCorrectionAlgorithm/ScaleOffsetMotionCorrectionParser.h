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
#ifndef _SCALEOFFSETMOTION_CORRECTION_PARSER_H_
#define _SCALEOFFSETMOTION_CORRECTION_PARSER_H_

#include <stdio.h>
#include <string>
#include "ScaleOffsetMotionStructures.h"



/**
 * @class ScaleOffsetCorrectionParser ScaleOffsetCorrectionParser.h ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionParser.h
 * @brief
 * @author
 * @date Nov 16, 2011
 * @version 1.0
 */
class ScaleOffsetCorrectionParser
{
  public:
    ScaleOffsetCorrectionParser();
    virtual ~ScaleOffsetCorrectionParser();

    int parseArguments(int argc, char **argv, TomoInputs* Input);
    int readParameterFile(const std::string &filepath,TomoInputs* ParsedInput,Sino* Sinogram,Geom* Geometry);
    void initializeSinoParameters(Sino* Sinogram,TomoInputs* ParsedInput);
    void initializeGeomParameters(Sino* Sinogram,Geom* Geometry,TomoInputs* ParsedInput);
    DATA_TYPE absMaxArray(DATA_TYPE* Array, uint16_t NumElts);

    /**
     * @brief Copys the std::string contents into a newly malloc'ed char array which
     * the programmer will need to free when they are finished with it.
     * @param fname The filename to copy
     */
    char* copyFilenameToNewCharBuffer( const std::string &fname);


  private:
    uint64_t startm;
    uint64_t stopm;

    ScaleOffsetCorrectionParser(const ScaleOffsetCorrectionParser&); // Copy Constructor Not Implemented
    void operator=(const ScaleOffsetCorrectionParser&); // Operator '=' Not Implemented

};

#endif //ComputationInputs
