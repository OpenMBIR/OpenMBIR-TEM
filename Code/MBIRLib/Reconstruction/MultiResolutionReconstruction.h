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


#ifndef MULTIRESOLUTIONSOC_H_
#define MULTIRESOLUTIONSOC_H_


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/FilterPipeline.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/Reconstruction/ReconstructionEngine.h"

/**
 * @brief This class controls the multiresolution reconstruction of an input
 * data set
 */
class MultiResolutionReconstruction : public FilterPipeline
{
  public:
    MXA_SHARED_POINTERS(MultiResolutionReconstruction)
    MXA_TYPE_MACRO_SUPER(MultiResolutionReconstruction, FilterPipeline)
    MXA_STATIC_NEW_MACRO(MultiResolutionReconstruction)

    virtual ~MultiResolutionReconstruction();

    /**
     * @brief Cancel the operation
     */
    virtual void setCancel(bool value);
    virtual bool getCancel();

    MXA_INSTANCE_PROPERTY(bool, Debug)

    MXA_INSTANCE_STRING_PROPERTY(InputFile)
    MXA_INSTANCE_STRING_PROPERTY(TempDir)
    MXA_INSTANCE_STRING_PROPERTY(OutputFile)
    MXA_INSTANCE_STRING_PROPERTY(BrightFieldFile)
    MXA_INSTANCE_STRING_PROPERTY(InitialReconstructionFile)
    MXA_INSTANCE_PROPERTY(bool, DeleteTempFiles)

    MXA_INSTANCE_PROPERTY(int, NumberResolutions)
    MXA_INSTANCE_PROPERTY(float, SampleThickness)
    MXA_INSTANCE_PROPERTY(float, TargetGain)
    MXA_INSTANCE_PROPERTY(float, BraggThreshold)
    MXA_INSTANCE_PROPERTY(float, BraggDelta)
    MXA_INSTANCE_PROPERTY(float, BfOffset)
    MXA_INSTANCE_PROPERTY(float, StopThreshold)
    MXA_INSTANCE_PROPERTY(int, OuterIterations)
    MXA_INSTANCE_PROPERTY(int, InnerIterations)
    MXA_INSTANCE_PROPERTY(float, SigmaX)
    MXA_INSTANCE_PROPERTY(float, MRFShapeParameter)
    MXA_INSTANCE_PROPERTY(float, DefaultOffsetValue)
    MXA_INSTANCE_PROPERTY(bool, UseDefaultOffset)
    MXA_INSTANCE_PROPERTY(int, FinalResolution)
    MXA_INSTANCE_PROPERTY(bool, ExtendObject)
    MXA_INSTANCE_PROPERTY(bool, InterpolateInitialReconstruction)
    MXA_INSTANCE_PROPERTY(float, DefaultVariance)
    MXA_INSTANCE_PROPERTY(float, InitialReconstructionValue)

    MXA_INSTANCE_PROPERTY(std::vector<float>, Tilts)
    MXA_INSTANCE_PROPERTY(AdvancedParametersPtr, AdvParams)


    /**
     * @brief If this vector is set, ie, length = 6, then we are going to use
     * a subvolume to reconstruct. The values are encoded in the vector as
     * xMin, yMin, zMin, xMax, yMax, zMax.
     */
    MXA_INSTANCE_PROPERTY(std::vector<uint16_t>, Subvolume)

    MXA_INSTANCE_PROPERTY(std::vector<uint8_t>, ViewMasks)


    /**
     * @brief
     */
    virtual void execute();


    /**
   * @brief Prints the values of the input variables
   * @param inputs
   * @param out
   */
    void printInputs(TomoInputsPtr inputs, std::ostream &out);

    /**
     * @brief Trys to estimate the amount of memory that will be used.
     * @param inputs
     * @param bf_inputs
     */
    void memCalculate(TomoInputsPtr inputs);//,TomoInputsPtr bf_inputs

  protected:
    MultiResolutionReconstruction();

  private:
    bool                 m_Cancel;
    ReconstructionEngine::Pointer   m_CurrentEngine;

    MultiResolutionReconstruction(const MultiResolutionReconstruction&); // Copy Constructor Not Implemented
    void operator=(const MultiResolutionReconstruction&); // Operator '=' Not Implemented
};

#endif /* MULTIRESOLUTIONSOC_H_ */
