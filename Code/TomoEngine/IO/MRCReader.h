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

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "MXA/Common/IO/MXAFileReader64.h"


#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/IO/MRCHeader.h"


/**
 * @class MRCReader MRCReader.h EIMTomo/IO/MRCReader.h
 * @brief A Class to read a .mrc file
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Nov 22, 2011
 * @version 1.0
 */
class TomoEngine_EXPORT MRCReader
{
  public:

    MXA_SHARED_POINTERS(MRCReader);
    static Pointer New(bool deleteMemory = true)
    {
      Pointer sharedPtr (new MRCReader(deleteMemory));
      return sharedPtr;
    }

    MXA_TYPE_MACRO(MRCReader)



    virtual ~MRCReader();

    MXA_INSTANCE_PROPERTY(bool, DeleteMemory);

    /**
     * @brief This method ONLY reads the header section of the file
     * @param filepath The path to the input file
     * @param header Pointer to an MRCHeader structure where the header values will
     * be stored.
     * @return Negative on Error.
     */
    int readHeader(const std::string &filepath, MRCHeader* header);

    /**
     * @brief Reads the entire file or a subpart of the file into memory. If a
     * subpart of the file is to be read then the voxelMin and voxelMax pointers
     * should be non-null and point to 3 element arrays describing the minimum and
     * maximum voxel indices along each axis that should be read. For example if
     * the dimensions of the volume are 10 x 20 x 30 (columns, rows, sections) and you want to read
     * a subvolume of (5,8) along x, (12, 15) along y and (23, 30) along z then
     * the programmer would pass in arrays of the following values:
     * @code
     * int min[3] = {5, 12, 23};
     * int max[3] = {8, 15, 30};
     * MRCReader reader;
     * reader.read(aFilePath, min, max);
     * @endcode
     *  Note that the ranges are INCLUSIVE of the max value, ie, if the number of columns
     *  has 20 voxels and you want to read all the voxels then the range is 0-19.
     *
     *  If the entire volume is requested then simply pass in NULL for the pointers.
     *
     * @param filepath The path to the input file.
     * @param voxelMin The minimum index value for the voxels [ x, y, z ]
     * @param voxelMax The maximun index value for the voxels [ x, y, z ]
     * @return Negative on Error.
     */
    int read(const std::string &filepath, int* voxelMin = NULL, int* voxelMax = NULL);

    bool readPartialVolume(MXAFileReader64 &reader, char* dataPtr,
                       size_t typeSize, size_t nVoxels,
                       int* voxelMin, int* voxelMax);

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

    std::string getLabelField(int index);

    void printHeader(MRCHeader* header, std::ostream &out);

  protected:
    MRCReader();
    /**
     * @brief Constructur
     * @param deleteMemory Should this class delete the memory allocated to hold
     * the voxel data.
     */
    explicit MRCReader(bool deleteMemory);

  private:
    MRCHeader* m_Header;

    uint8_t*   m_UInt8Data;
    int16_t*   m_Int16Data;
    uint16_t*  m_UInt16Data;
    float*     m_FloatData;



    MRCReader(const MRCReader&); // Copy Constructor Not Implemented
    void operator=(const MRCReader&); // Operator '=' Not Implemented
};

#endif /* MRCREADER_H_ */
