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

#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"

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



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectorResponse::execute()
{


  DATA_TYPE r,sum=0,rmin,ProfileCenterR,ProfileCenterT,TempConst,tmin;
  DATA_TYPE r0 = -(m_BeamWidth)/2;
  DATA_TYPE StepSize = m_BeamWidth/BEAM_RESOLUTION;
  int16_t i,j,k,p,ProfileIndex;
  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr inputs = getTomoInputs();

  size_t dims[3] = {1, sinogram->N_theta,DETECTOR_RESPONSE_BINS};
  RealVolumeType::Pointer H = RealVolumeType::New(dims, "DetectorResponse");

  //H = (DATA_TYPE***)get_3D(1, m_Sinogram->N_theta,DETECTOR_RESPONSE_BINS, sizeof(DATA_TYPE));//change from 1 to DETECTOR_RESPONSE_BINS
  TempConst=(PROFILE_RESOLUTION)/(2*inputs->delta_xz);

  for(k = 0 ; k < sinogram->N_theta; k++)
  {
    for (i = 0; i < DETECTOR_RESPONSE_BINS; i++) //displacement along r
    {
      ProfileCenterR = i*m_OffsetR;
      rmin = ProfileCenterR - inputs->delta_xz;
      for (j = 0 ; j < 1; j++)//displacement along t ;change to DETECTOR_RESPONSE_BINS later
      {
        ProfileCenterT = j*m_OffsetT;
        tmin = ProfileCenterT - inputs->delta_xy/2;
        sum = 0;
        for (p=0; p < BEAM_RESOLUTION; p++)
        {
          r = r0 + p*StepSize;
          if(r < rmin)
            continue;

          ProfileIndex = (int32_t)floor((r - rmin) * TempConst);
          if(ProfileIndex < 0)
          {
            ProfileIndex = 0;
          }
          if(ProfileIndex >= PROFILE_RESOLUTION)
          {
            ProfileIndex = PROFILE_RESOLUTION - 1;
          }
          //sum += (m_VoxelProfile->d[k][ProfileIndex] * m_BeamProfile->d[p]);
          sum += (m_VoxelProfile->getValue(k,ProfileIndex) * m_BeamProfile->d[p]);
        }
        H->d[j][k][i] = sum;
      }
    }
  }

  // Set a reference to the data in this class so the calling class can get the response
  setResponse(H);

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Calculating the Detector Response", 0, UpdateProgressMessage);
}
