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

#include "DetectorResponse.h"

#include "MBIRLib/Common/EIMMath.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectorResponse::~DetectorResponse()
{
  m_Response = RealVolumeType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectorResponse::DetectorResponse()
{
}


#define MAKE_2D_INDEX(index, dim1, idx0, idx1)\
  index = ((dim1) * (idx0)) + (idx1);

#define MAKE_3D_INDEX(index, dim1, dim2, idx0, idx1, idx2) \
  index = ((dim1)*(dim2)*(idx0)) + ((dim2)*(idx1)) +(idx2);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectorResponse::execute()
{
  int16_t i, j, k, ProfileIndex;
  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr inputs = getTomoInputs();
  AdvancedParametersPtr advParams = getAdvParams();

  Real_t beamWidth = m_DetectorParameters->getBeamWidth();
  Real_t offsetR = m_DetectorParameters->getOffsetR();
  Real_t offsetT = m_DetectorParameters->getOffsetT();
  RealArrayType::Pointer beamProfile = m_DetectorParameters->getBeamProfile();


  Real_t r, sum = 0, rmin, ProfileCenterR, ProfileCenterT, TempConst, tmin;
  Real_t r0 = -(beamWidth) / 2;
  Real_t StepSize = beamWidth / advParams->BEAM_RESOLUTION;
  size_t dims[3] = {1, sinogram->N_theta, advParams->DETECTOR_RESPONSE_BINS};
  RealVolumeType::Pointer H = RealVolumeType::New(dims, "DetectorResponse");

  //H = (DATA_TYPE***)get_3D(1, m_Sinogram->N_theta,DETECTOR_RESPONSE_BINS, sizeof(DATA_TYPE));//change from 1 to DETECTOR_RESPONSE_BINS
  TempConst = (advParams->PROFILE_RESOLUTION) / (2 * inputs->delta_xz);

  for(k = 0 ; k < sinogram->N_theta; k++)
  {
    for (i = 0; i < advParams->DETECTOR_RESPONSE_BINS; i++) //displacement along r
    {
      ProfileCenterR = i * offsetR;
      rmin = ProfileCenterR - inputs->delta_xz;
      for (j = 0; j < 1; j++) //displacement along t ;change to DETECTOR_RESPONSE_BINS later
      {
        ProfileCenterT = j * offsetT;
        tmin = ProfileCenterT - inputs->delta_xy / 2;
        sum = 0;
        for (uint32_t p = 0; p < advParams->BEAM_RESOLUTION; p++)
        {
          r = r0 + p * StepSize;
          if(r < rmin) { continue; }

          ProfileIndex = (int32_t)floor((r - rmin) * TempConst);
          if(ProfileIndex < 0)
          {
            ProfileIndex = 0;
          }
          if(ProfileIndex >= advParams->PROFILE_RESOLUTION)
          {
            ProfileIndex = advParams->PROFILE_RESOLUTION - 1;
          }
          sum += m_VoxelProfile->getValue(k, ProfileIndex) * beamProfile->d[p];
//          sum += (m_VoxelProfile->d[k][ProfileIndex] * m_BeamProfile->d[p]);//;*BeamProfile[l]);
        }
        H->setValue(sum, j, k, i);
      }
    }
  }

  // Set a reference to the data in this class so the calling class can get the response
  setResponse(H);

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Calculating the Detector Response", 0, UpdateProgressMessage);
}
