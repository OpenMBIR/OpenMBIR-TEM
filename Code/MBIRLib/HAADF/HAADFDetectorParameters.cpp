/*
 * HAADFParameters.cpp
 *
 *  Created on: Sep 5, 2012
 *      Author: mjackson
 */

#include "HAADFDetectorParameters.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFDetectorParameters::HAADFDetectorParameters()
{
  // TODO Auto-generated constructor stub

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFDetectorParameters::~HAADFDetectorParameters()
{
  // TODO Auto-generated destructor stub
}



/* Initializes the global variables cosine and sine to speed up computation
 */
void HAADFDetectorParameters::calculateSinCos(SinogramPtr m_Sinogram)
{
  uint16_t i;
  size_t dims[1] = { m_Sinogram->N_theta };
  m_cosine = RealArrayType::New(dims, "cosine");
  m_sine = RealArrayType::New(dims, "sine");

  for(i=0;i<m_Sinogram->N_theta;i++)
  {
    m_cosine->d[i]=cos(m_Sinogram->angles[i]);
    m_sine->d[i]=sin(m_Sinogram->angles[i]);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFDetectorParameters::initializeBeamProfile(SinogramPtr m_Sinogram, AdvancedParametersPtr m_AdvParams)
{
  uint16_t i;
  Real_t sum=0,W;
  size_t dims[1] = {m_AdvParams->BEAM_RESOLUTION };
  m_BeamProfile = RealArrayType::New(dims, "BeamProfile");
  W=m_BEAM_WIDTH/2;
  for (i=0; i < m_AdvParams->BEAM_RESOLUTION ;i++)
  {
    m_BeamProfile->d[i] = 0.54 - 0.46*cos((2.0*M_PI/m_AdvParams->BEAM_RESOLUTION)*i);
    sum=sum+m_BeamProfile->d[i];
  }

  //Normalize the beam to have an area of 1
  for (i=0; i <m_AdvParams->BEAM_RESOLUTION ;i++)
  {
    m_BeamProfile->d[i]/=sum;
    m_BeamProfile->d[i]/=m_Sinogram->delta_t;//This is for proper normalization
    // printf("%lf\n",BeamProfile->d[i]);
  }



}
