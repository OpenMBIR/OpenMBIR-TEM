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
#ifndef RAWGEOMETRYWRITER_H_
#define RAWGEOMETRYWRITER_H_


#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"

/**
 * @class
 * @brief
 * @author
 * @date
 * @version
 */
class TomoEngine_EXPORT RawGeometryWriter : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(RawGeometryWriter)
    MXA_STATIC_NEW_MACRO(RawGeometryWriter)
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, RawGeometryWriter)
    MXA_TYPE_MACRO_SUPER(RawGeometryWriter, TomoFilter)

    virtual ~RawGeometryWriter();

    MXA_INSTANCE_STRING_PROPERTY(FilePath)

    void execute();


  protected:
    RawGeometryWriter();

  private:
    RawGeometryWriter(const RawGeometryWriter&); // Copy Constructor Not Implemented
    void operator=(const RawGeometryWriter&); // Operator '=' Not Implemented
};

#endif /* RAWGEOMETRYWRITER_H_ */
