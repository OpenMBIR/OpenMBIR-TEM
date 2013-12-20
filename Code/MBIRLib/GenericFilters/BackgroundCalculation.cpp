/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
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

#include "BackgroundCalculation.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BackgroundCalculation::BackgroundCalculation(std::string filePath, int x, int y, int width, int height) :
m_Path(filePath)
{
    m_Reader = MRCReader::New();
    MRCHeader* header = new MRCHeader();
    
    
    if ( !m_Reader->readHeader(m_Path, header) )
    {
        std::cout << "m_Reader did not read header from filepath correctly!\n";
        return;
    }
    
    int z = header->nz;
    int mode = header->mode;
    
    int min[3] = {x, y, 0};
    int max[3] = {x+width, y+height, z-1};
    int nVoxels = width*height*z;
    
    if ( !m_Reader->read(m_Path, min, max) )
    {
        std::cout << "m_Reader did not read filepath correctly!\n";
        return;
    }
    
    switch (mode)
    {
        case TYPE_BYTES:
        {
            if (header->imodFlags == 1)
            {
                // Signed bytes
                computeMean<signed char>(m_Reader->getDataPointer(), nVoxels, TYPE_SIGNED_BYTES);
            }
            else
            {
                // Unsigned bytes
                computeMean<unsigned char>(m_Reader->getDataPointer(), nVoxels, TYPE_UNSIGNED_BYTES);
            }
            break;
        }
        case TYPE_SIGNED_SHORT_INT:
        {
            computeMean<signed short int>(m_Reader->getDataPointer(), nVoxels, TYPE_SIGNED_SHORT_INT);
            break;
        }
        case TYPE_FLOAT:
        {
            computeMean<float>(m_Reader->getDataPointer(), nVoxels, TYPE_FLOAT);
            break;
        }
    }

    delete header;
    header = NULL;
    
    /* Types of pixel in image.  Values used by IMOD:
        0 = unsigned or signed bytes depending on flag in imodStamp
            only unsigned bytes before IMOD 4.2.23 (iModFlags == 1, signed)
        1 = signed short integers (16 bits)
        2 = float
        3 = short * 2, (used for complex data)
        4 = float * 2, (used for complex data)
        6 = unsigned 16-bit integers (non-standard)
        16 = unsigned char * 3 (for rgb data, non-standard) */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BackgroundCalculation::~BackgroundCalculation()
{
    
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template<class T> void BackgroundCalculation::computeMean(void *ptr, int numVoxels, Type type)
{
    double mean = 0.0;
    double sum = 0.0;
    
    if (type == TYPE_SIGNED_BYTES)
    {
        signed char* dataPtr = (signed char*)ptr;
        for (int i=0; i<numVoxels; i++)
        {
            sum = sum + dataPtr[i];
        }
    }
    else if (type == TYPE_UNSIGNED_BYTES)
    {
        unsigned char* dataPtr = (unsigned char*)ptr;
        for (int i=0; i<numVoxels; i++)
        {
            sum = sum + dataPtr[i];
        }
    }
    else if (type == TYPE_SIGNED_SHORT_INT)
    {
        signed short int* dataPtr = (signed short int*)ptr;
        for (int i=0; i<numVoxels; i++)
        {
            sum = sum + dataPtr[i];
        }
    }
    else if (type == TYPE_FLOAT)
    {
        float* dataPtr = (float*)ptr;
        for (int i=0; i<numVoxels; i++)
        {
            sum = sum + dataPtr[i];
        }
    }
    else
    {
        std::cout << "BackgroundCalculation::computeMean(...) - No type recognized\n";
        m_Mean = 0.0;
        return;
    }
    
    mean = sum/numVoxels;
    
    m_Mean = mean;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double BackgroundCalculation::getMeanValue()
{
    return m_Mean;
}



