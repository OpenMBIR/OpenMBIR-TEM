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

#include "HAADFDetectorParameters.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFDetectorParameters::HAADFDetectorParameters()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFDetectorParameters::~HAADFDetectorParameters()
{
}

// -----------------------------------------------------------------------------
/* Initializes the global variables cosine and sine to speed up computation
 */
void HAADFDetectorParameters::calculateSinCos(SinogramPtr m_Sinogram)
{
  uint16_t i;
  size_t dims[1] =
  { m_Sinogram->N_theta };
  m_cosine = RealArrayType::New(dims, "cosine");
  m_sine = RealArrayType::New(dims, "sine");

  for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    m_cosine->d[i] = cos(m_Sinogram->angles[i]);
    m_sine->d[i] = sin(m_Sinogram->angles[i]);
  }
}

// -----------------------------------------------------------------------------
// Beam profile used for A-matrix
// -----------------------------------------------------------------------------
void HAADFDetectorParameters::initializeBeamProfile(SinogramPtr m_Sinogram, AdvancedParametersPtr m_AdvParams)
{
  uint16_t i;
  Real_t sum = 0, W;
  size_t dims[1] =
  { m_AdvParams->BEAM_RESOLUTION };
  m_BeamProfile = RealArrayType::New(dims, "BeamProfile");
  W = m_BeamWidth / 2;
  for (i = 0; i < m_AdvParams->BEAM_RESOLUTION; i++)
  {
    m_BeamProfile->d[i] = 0.54 - 0.46 * cos((2.0 * M_PI / m_AdvParams->BEAM_RESOLUTION) * i);
    sum = sum + m_BeamProfile->d[i];
  }

  //Normalize the beam to have an area of 1
  for (i = 0; i < m_AdvParams->BEAM_RESOLUTION; i++)
  {
    m_BeamProfile->d[i] /= sum;
    m_BeamProfile->d[i] /= m_Sinogram->delta_t; //This is for proper normalization
    // printf("%lf\n",BeamProfile->d[i]);
  }

}
