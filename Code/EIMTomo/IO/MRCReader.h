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

#ifndef MRCREADER_H_
#define MRCREADER_H_

#include <string>

#include "MRCHeader.h"



/**
 * @class MRCReader MRCReader.h EIMTomo/IO/MRCReader.h
 * @brief A Class to read a .mrc file
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Nov 22, 2011
 * @version 1.0
 */
class MRCReader
{
  public:
    /**
     * @brief Constructur
     * @param deleteMemory Should this class delete the memory allocated to hold
     * the voxel data.
     */
    explicit MRCReader(bool deleteMemory = true);
    virtual ~MRCReader();

    /**
     * @brief This method ONLY reads the header section of the file
     * @param filepath The path to the input file
     * @param header Pointer to an MRCHeader structure where the header values will
     * be stored.
     * @return Negative on Error.
     */
    int readHeader(const std::string &filepath, MRCHeader* header);

    /**
     * @brief Reads the entire file into memory
     * @param filepath The path to the input file.
     * @return Negative on Error.
     */
    int read(const std::string &filepath);

    /**
     * @brief Returns the pointer to the data from the file.
     */
    void* getDataPointer();

    /**
     * @brief Returns the header structure. Note that this pointer is owned by
     * this class and will be deleted when this class is destroyed.
     * @return
     */
    MRCHeader* getHeader();

  protected:
    MRCReader();

  private:
    MRCHeader* m_Header;

    uint8_t*  m_UInt8Data;
    int16_t*  m_Int16Data;
    uint16_t*  m_UInt16Data;
    float*    m_FloatData;
    bool      m_DeleteMemory;

    MRCReader(const MRCReader&); // Copy Constructor Not Implemented
    void operator=(const MRCReader&); // Operator '=' Not Implemented
};

#endif /* MRCREADER_H_ */
