/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
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
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
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

#ifndef MRCWRITER_H_
#define MRCWRITER_H_


#include <string>

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/GenericFilters/TomoFilter.h"
#include "MBIRLib/IOFilters/MRCHeader.h"

/**
 * @class MRCReader MRCReader.h EIMTomo/IO/MRCWriter.h
 * @brief A Class to write a .mrc file
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Feb 22, 2012
 * @version 1.0
 */
class MBIRLib_EXPORT MRCWriter : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(MRCWriter)
    static Pointer New()
    {
      Pointer sharedPtr (new MRCWriter);
      return sharedPtr;
    }

    MXA_TYPE_MACRO(MRCWriter)

    virtual ~MRCWriter();

    MXA_INSTANCE_PROPERTY(bool, DeleteMemory)
    MXA_INSTANCE_STRING_PROPERTY(OutputFile)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)
    MXA_INSTANCE_PROPERTY(MRCHeader*, MRCHeader)
    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, XDims)
    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, YDims)
    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, ZDims)

    void execute();

    void initializeMRCHeader(MRCHeader* header);

    /**
     * @brief This method ONLY reads the header section of the file
     * @param filepath The path to the input file
     * @param header Pointer to an MRCHeader structure where the header values will
     * be stored.
     * @return Negative on Error.
     */
    int writeHeader();

    int write();

  protected:

    /**
     * @brief Constructur
     */
    explicit MRCWriter();

  private:
    MRCWriter(const MRCWriter&); // Copy Constructor Not Implemented
    void operator=(const MRCWriter&); // Operator '=' Not Implemented

};

#endif /* MRCWRITER_H_ */
