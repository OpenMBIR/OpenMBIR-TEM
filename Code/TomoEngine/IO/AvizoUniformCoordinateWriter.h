/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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

#ifndef AvizoUniformCoordinateWriter_H_
#define AvizoUniformCoordinateWriter_H_

#include <string>

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"
#include "MXA/Common/IO/MXAFileWriter64.h"


/**
 * @class AvizoUniformCoordinateWriter AvizoUniformCoordinateWriter.h DREAM3DLib/IOFilters/AvizoUniformCoordinateWriter.h
 * @brief Writes out a native Avizo Uniform Coordinate file
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Aug 9, 2012
 * @version 1.0
 */
class AvizoUniformCoordinateWriter : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(AvizoUniformCoordinateWriter);
    MXA_STATIC_NEW_MACRO(AvizoUniformCoordinateWriter);
    MXA_TYPE_MACRO_SUPER(AvizoUniformCoordinateWriter, AbstractFilter);

    virtual ~AvizoUniformCoordinateWriter();

    MXA_INSTANCE_PROPERTY(bool, DeleteMemory)
    MXA_INSTANCE_STRING_PROPERTY(OutputFile)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)
    MXA_INSTANCE_PROPERTY(bool, WriteBinaryFile)

    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, XDims)
    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, YDims)
    MXA_INSTANCE_VEC2_PROPERTY(uint16_t, ZDims)


    void execute();

    int write();

  protected:
    AvizoUniformCoordinateWriter();

    /**
     * @brief Generates the Avizo Header for this file
     * @return The header as a string
     */
    std::string generateHeader();

    /**
     * @brief Writes the data to the Avizo file
     * @param writer The MXAFileWriter object
     * @return Error code
     */
    int writeData(MXAFileWriter64 &writer);

  private:
    int32_t* m_GrainIds;

    AvizoUniformCoordinateWriter(const AvizoUniformCoordinateWriter&); // Copy Constructor Not Implemented
    void operator=(const AvizoUniformCoordinateWriter&); // Operator '=' Not Implemented
};


#endif
