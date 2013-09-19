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
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Purdue
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



#ifndef SCALEOFFSETCORRECTIONCONSTANTS_H_
#define SCALEOFFSETCORRECTIONCONSTANTS_H_

#include <string>


#define INDEX_3(i, j, k)\
    ((9*(i)) + (3*(j)) + ((k)))

namespace VoxelUpdateType
{
  const unsigned int RegularRandomOrderUpdate = 0;
  const unsigned int HomogeniousUpdate = 1;
  const unsigned int NonHomogeniousUpdate = 2;
}



//#define FORWARD_PROJECT_MODE


//Region Of Interest for calculating the stopping criteria. Should be on with stopping threshold
#define ROI 1
#define EIMTOMO_USE_QGGMRF 1 //Was used earlier to switch priors
#define BRAGG_CORRECTION //TODO : Get rid and make this the default case 

//#define WRITE_INTERMEDIATE_RESULTS
//#define COST_CALCULATE
//#define DEBUG_COSTS
#define DEBUG
//#define RANDOM_ORDER_UPDATES

#define DefBraggThreshold 5e10//Thresold in normalized uints. This value ensures at the start we do a "regular" BF recon
#define BF_OFFSET 32768//23696 //- Bio data set
#define BF_MAX  1//1865//30369//5689//42122//1865//32768//
#define BF_T 3e10 //Ratio of model mismatch to noise variance
//#define IDENTITY_NOISE_MODEL
#define POSITIVITY_CONSTRAINT
#define NHICD
#define NUM_NON_HOMOGENOUS_ITER 20
#define SUB_ITER 1
#define NUM_HOM_ITER 20

namespace ScaleOffsetCorrection
{

//  const std::string OutputDirectory("ScaleOffsetCorrection_Output");
  const std::string DetectorResponseFile("DetectorResponse.bin");
  const std::string CostFunctionFile("CostFunc.bin");
  const std::string FinalGainParametersFile("FinalGainParameters.bin");
  const std::string FinalOffsetParametersFile("FinalOffsetParameters.bin");
  const std::string FinalVariancesFile("FinalVariances.bin");

  const std::string VoxelProfileFile("VoxelProfile.bin");
  const std::string ForwardProjectedObjectFile("ForwardProjectedObject.bin");
  const std::string CostFunctionCoefficientsFile("CostFunctionCoefficients.bin");
  const std::string FilteredMagMapFile("FilteredMagMap.bin");
  const std::string MagnitudeMapFile("MagnitudeMap.bin");

  const std::string ReconstructedVtkFile("ReconstructedVolume.vtk");
  const std::string ReconstructedSinogramFile("ReconstructedSinogram.bin");
  const std::string ReconstructedObjectFile("ReconstructedObject.bin");
  const std::string ReconstructedMrcFile("ReconstructedVolume.rec");

  const std::string UpsampledBinFile("UpsampledObject.bin");

  namespace VTK
  {
    const std::string TomoVoxelScalarName("TomoVoxel");
  }

}



#endif /* SCALEOFFSETCORRECTIONCONSTANTS_H_ */
