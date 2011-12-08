/*
 * DetectorResponse.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: mjackson
 */

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
  Sinogram* sinogram = getSinogram();
  TomoInputs* inputs = getInputs();

  size_t dims[3] = {1, sinogram->N_theta,DETECTOR_RESPONSE_BINS};
  RealVolumeType::Pointer H = RealVolumeType::New(dims);
  H->setName("DetectorResponse");

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
          sum += (m_VoxelProfile->d[k][ProfileIndex] * m_BeamProfile->d[p]);//;*BeamProfile[l]);
        }
        H->d[j][k][i] = sum;
      }
    }
  }

  // Set a reference to the data in this class so the calling class can get the response
  setResponse(H);

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Done Calculating the Detector Response", 0, UpdateProgressMessage);
}
