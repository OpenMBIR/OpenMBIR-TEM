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
//#define FORWARD_PROJECT_MODE
#define EXTEND_OBJECT
#define X_SHRINK_FACTOR 0.6
//#define TARGET_GAIN 7872950.013398

#define X_STRETCH 1
#define Z_STRETCH 2

#define INTERPOLATE_FACTOR 24
//#define DEBUG ,
#define PROFILE_RESOLUTION 1536
//#define PI 4*atan(1)//3.14159265
//Beam Parameters - This is set to some number <<< Sinogram->delta_r.
//#define BEAM_WIDTH 0.050000
#define BEAM_RESOLUTION 512
#define AREA_WEIGHTED

//Region Of Interest for calculating the stopping criteria. Should be on with stopping threshold
#define ROI 1

#define STOPPING_THRESHOLD 0.001
#define STOPPING_THRESHOLD_I_k 0.001
#define STOPPING_THRESHOLD_Delta_k 0.001
#define STOPPING_THRESHOLD_Var_k 0.001
#define THRESHOLD_REDUCTION_FACTOR 1 //Dynamically lower the threshold by this amount. Set to 1 for no reduction
#define SURROGATE_FUNCTION
#define QGGMRF
#define QGGMRF_ITER 1
//#define DISTANCE_DRIVEN
//#define CORRECTION
//#define WRITE_INTERMEDIATE_RESULTS
#define COST_CALCULATE
#define DETECTOR_RESPONSE_BINS 64
#define JOINT_ESTIMATION
//#define GEOMETRIC_MEAN_CONSTRAINT
#define NOISE_MODEL
#define POSITIVITY_CONSTRAINT
//#define CIRCULAR_BOUNDARY_CONDITION
//#define DEBUG_CONSTRAINT_OPT
#define RANDOM_ORDER_UPDATES
//#define BRIGHT_FIELD
//#define NHICD 
#define NUM_NON_HOMOGENOUS_ITER 20

//#define FORWARD_PROJECT_MODE //this Flag just takes the input file , forward projects it and exits


namespace ScaleOffsetCorrection
{

//  const std::string OutputDirectory("ScaleOffsetCorrection_Output");
  const std::string DetectorResponseFile("DetectorResponse.bin");
  const std::string CostFunctionFile("CostFunc.bin");
  const std::string FinalGainParametersFile("FinalGainParameters.bin");
  const std::string FinalOffsetParametersFile("FinalOffsetParameters.bin");
  const std::string FinalVariancesFile("FinalVariances.bin");
  const std::string ReconstructedSinogramFile("ReconstructedSino.bin");
  const std::string VoxelProfileFile("VoxelProfile.bin");
  const std::string ForwardProjectedObjectFile("ForwardProjectedObject.bin");
  const std::string CostFunctionCoefficientsFile("CostFunctionCoefficients.bin");
  const std::string FilteredMagMapFile("FilteredMagMap.bin");
  const std::string MagnitudeMapFile("MagnitudeMap.bin");
  const std::string VtkGeometryFile("ReconstructedData.vtk");

  namespace VTK
  {
    const std::string TomoVoxelScalarName("TomoVoxel");
  }

}



#endif /* SCALEOFFSETCORRECTIONCONSTANTS_H_ */
