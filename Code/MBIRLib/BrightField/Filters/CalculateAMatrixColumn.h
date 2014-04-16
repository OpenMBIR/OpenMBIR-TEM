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
#ifndef CALCULATEAMATRIX_H_
#define CALCULATEAMATRIX_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/GenericFilters/TomoFilter.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/BrightField/BFAMatrixCol.h"

/*
 *
 */
class MBIRLib_EXPORT CalculateAMatrixColumn : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(CalculateAMatrixColumn)
    MXA_STATIC_NEW_MACRO(CalculateAMatrixColumn);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, CalculateAMatrixColumn);
    MXA_TYPE_MACRO_SUPER(CalculateAMatrixColumn, TomoFilter)

    virtual ~CalculateAMatrixColumn();


    MXA_INSTANCE_PROPERTY_OLD(Real_t*, Cosine, cosine);
    MXA_INSTANCE_PROPERTY_OLD(Real_t*, Sine, sine);

    MXA_INSTANCE_PROPERTY_OLD(Real_t, BeamWidth, BEAM_WIDTH);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, row, row);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, col, col);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, Slice, slice);
    MXA_INSTANCE_PROPERTY_OLD(Real_t**, VoxelProfile, VoxelProfile);
    MXA_INSTANCE_PROPERTY_OLD(Real_t*, D1, d1);
    MXA_INSTANCE_PROPERTY_OLD(Real_t*, D2, d2);
    MXA_INSTANCE_PROPERTY_OLD(BFAMatrixCol::Pointer, AMatrixCol, Ai)


    virtual void execute();

  protected:
    CalculateAMatrixColumn();
};

#endif /* CALCULATEAMATRIX_H_ */
