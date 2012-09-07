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


#ifndef _HAADFAMATRIXCOL_H_
#define _HAADFAMATRIXCOL_H_

#include "MXA/Common/MXASetGetMacros.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/HAADF/HAADFDetectorParameters.h"



/**
 * @brief Holds a column specific to the HAADF algorihms and inputs codes of the
 * A Matrix for Tomographic reconstructions
 */
class HAADFAMatrixCol
{
  public:

    MXA_SHARED_POINTERS(HAADFAMatrixCol);
    static Pointer New(size_t* dims, int32_t count)
    {
      Pointer sharedPtr (new HAADFAMatrixCol(dims, count));
      return sharedPtr;
    }

    virtual ~HAADFAMatrixCol();

    RealArrayType::Pointer valuesPtr;
    Real_t* values;
    UInt32ArrayType::Pointer indexPtr;
    uint32_t* index;
    uint64_t d0;
    uint64_t d1;

    uint32_t count; //The number of non zero values present in the column
    void setCount(uint32_t c);

    static HAADFAMatrixCol::Pointer calculateHAADFAMatrixColumnPartial(SinogramPtr sinogram,
                                                                       GeometryPtr geometry,
                                                                       TomoInputsPtr tomoInputs,
                                                                       AdvancedParametersPtr advParams,
                                                                       uint16_t row,
                                                                       uint16_t col,
                                                                       uint16_t slice,
                                                                       RealVolumeType::Pointer detectorResponse,
                                                                       HAADFDetectorParameters::Pointer haadfParameters);

  protected:
    HAADFAMatrixCol(size_t* dims, int32_t c);


  private:

    HAADFAMatrixCol(const HAADFAMatrixCol&); // Copy Constructor Not Implemented
    void operator=(const HAADFAMatrixCol&); // Operator '=' Not Implemented
};
#endif /* _HAADFAMATRIXCOL_H_ */
