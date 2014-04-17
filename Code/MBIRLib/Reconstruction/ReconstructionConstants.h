/* ============================================================================
 * Copyright (c) 2014 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2014 Singanallur Venkatakrishnan (Purdue University)
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



#ifndef _ReconstructionConstants_H_
#define _ReconstructionConstants_H_

#include <string>

#define INDEX_3(i, j, k)\
  ((9*(i)) + (3*(j)) + ((k)))


typedef double Real_t;

namespace SOC {
  enum TiltSelection {
    A_Tilt,
    B_Tilt
  };
}


namespace MBIR
{



  namespace VoxelUpdateType
  {
    const unsigned int RegularRandomOrderUpdate = 0;
    const unsigned int HomogeniousUpdate = 1;
    const unsigned int NonHomogeniousUpdate = 2;
  }



  namespace Constants
  {
    const unsigned int k_NumHomogeniousIter = 20;
    const unsigned int k_NumNonHomogeniousIter = 20;
    const Real_t k_QGGMRF_Gamma = 5.0;
    const Real_t k_MaxAngleStretch = 75.0;
  }


  namespace Defaults
  {

    //  const std::string OutputDirectory("ScaleOffsetCorrection_Output");
    const std::string DetectorResponseFile("DetectorResponse.bin");
    const std::string CostFunctionFile("CostFunc.bin");
    const std::string FinalGainParametersFile("FinalGainParameters.bin");
    const std::string FinalOffsetParametersFile("FinalOffsetParameters.bin");
    const std::string FinalVariancesFile("FinalVariances.bin");

    const std::string VoxelProfileFile("VoxelProfile.bin");
    const std::string BFHAADF_ForwardProjectedObjectFile("BFHAADF_ForwardProjectedObject.bin");
    const std::string CostFunctionCoefficientsFile("CostFunctionCoefficients.bin");
    const std::string FilteredMagMapFile("FilteredMagMap.bin");
    const std::string MagnitudeMapFile("MagnitudeMap.bin");
    const std::string BraggSelectorFile("BraggSelector.mrc");

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

}


#endif /* _ReconstructionConstants_H_ */
