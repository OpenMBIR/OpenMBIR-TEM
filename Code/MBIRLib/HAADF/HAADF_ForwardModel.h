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

/*
 * The gains,offsets, variance (the very existence of those) + calculation of the
 * AMatrixColumn + Voxel Line Response ~= Forward Model specification
 *
 *
 * So are the Weights. The way we set the weight is a function of the forward model
 * too. And the nuisance parameters are also part of the forward model.
 *
 * The forward model is used in a bunch of ways. Initially it forward projects the initial
 * volume and computes a error between the data and the "forward projection".

 * Also in the iterations it is used to compute the various parameters (THETA1 and THETA2)
 * both of which require scanning through a portion of the the AMatrix, Voxel Line response,
 * Nuisance parameters , weights, bright field data etc. Which reminds me the BrightField
 * (which we have disable in the GUI is also a part of the forward model if its available)
 *
 * Not to the AMatrix and Voxel Line but to the nuisance parameters. So in some way
 * our forward model has some "unknown" things hence forcing us to update some values
 * of it. It like having an equation some coefficients I know , some I am going to figure out.
 *
 * So "Forward Project" is the algorithm and the "Forward Model" is the input into
 * the algorithm and the output from the algorithm is a reconstructed volume that
 * we compare to the original MRC data which gives us the "Error Sinogram"?
 *
 * Right. The output is not a volume but a tilt series (so we can compare it to
 * the acquired data).
 */
#ifndef _HAADF_ForwardModel_H_
#define _HAADF_ForwardModel_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/AbstractFilter.h"
#include "MBIRLib/Common/Observer.h"
#include "MBIRLib/Common/TomoArray.hpp"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/GenericFilters/CostData.h"
#include "MBIRLib/GenericFilters/DetectorParameters.h"

class HAADF_ForwardModel : public Observable
{

  public:
    MXA_SHARED_POINTERS(HAADF_ForwardModel)
    MXA_TYPE_MACRO(HAADF_ForwardModel)
    MXA_STATIC_NEW_MACRO(HAADF_ForwardModel)

    virtual ~HAADF_ForwardModel();

    MXA_INSTANCE_PROPERTY(bool, Verbose)
    MXA_INSTANCE_PROPERTY(bool, VeryVerbose)
    MXA_INSTANCE_PROPERTY(int, ErrorCondition)
    MXA_INSTANCE_PROPERTY(bool, Cancel)


    MXA_INSTANCE_PROPERTY(AdvancedParametersPtr, AdvParams)
    MXA_INSTANCE_PROPERTY(TomoInputsPtr, TomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, Sinogram)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)

    // Bright Field Data Flag
    MXA_INSTANCE_PROPERTY(Real_t, TargetGain)
    MXA_INSTANCE_PROPERTY(bool, BF_Flag)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialGain)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialOffset)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialVariance)

    // These are the Nuisance Parameters that we need to solve for
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, I_0) //Gains
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Mu) //Offset
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Alpha) //Noise variance refinement factor


    void writeNuisanceParameters(SinogramPtr sinogram);
    void writeSinogramFile(SinogramPtr sinogram,
                           RealVolumeType::Pointer Final_Sinogram);
    void writeReconstructionFile(const std::string &filepath);
    void writeVtkFile(const std::string &vtkFile, uint16_t cropStart, uint16_t cropEnd);
    void writeMRCFile(const std::string &mrcFile, uint16_t cropStart, uint16_t cropEnd);
    void writeAvizoFile(const std::string &file, uint16_t cropStart, uint16_t cropEnd);

    int createInitialGainsData();
    int createInitialOffsetsData();
    int createInitialVariancesData();

    void initializeHt(RealVolumeType::Pointer H_t, Real_t OffsetT);


  protected:
    HAADF_ForwardModel();

  private:


    HAADF_ForwardModel(const HAADF_ForwardModel&); // Copy Constructor Not Implemented
    void operator=(const HAADF_ForwardModel&); // Operator '=' Not Implemented
};

#endif /* _HAADF_ForwardModel_H_ */
