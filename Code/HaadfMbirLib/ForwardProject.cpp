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


#include "ForwardProject.h"

#include <sstream>


#include "ReconstructionCoreLib/Common/EIMMath.h"
#include "HaadfMbirLib/HaadfMbirConstants.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ForwardProject::ForwardProject(Sinogram* sinogram,
                               Geometry* geometry,
                               std::vector<AMatrixCol::Pointer> &tempCol,
                               std::vector<AMatrixCol::Pointer> &voxelLineResponse,
                               RealVolumeType::Pointer yEst,
                               ScaleOffsetParams* nuisanceParams,
                               uint16_t tilt,
                               Observable* obs) :
                  m_Sinogram(sinogram),
                  m_Geometry(geometry),
                  TempCol(tempCol),
                  VoxelLineResponse(voxelLineResponse),
                  Y_Est(yEst),
                  NuisanceParams(nuisanceParams),
                  m_Tilt(tilt),
                  m_Observable(obs)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ForwardProject::~ForwardProject()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ForwardProject::operator()() const
{
  std::stringstream ss;
  ss << "Forward projecting Z-Slice " << m_Tilt << "/" << m_Geometry->N_z;
  if (NULL != m_Observable)
  {
    m_Observable->notify(ss.str(), 0, Observable::UpdateProgressMessage);
  }

  uint32_t j = m_Tilt;
  uint16_t VoxelLineAccessCounter;
  uint32_t Index;
  Real_t ttmp = 0.0;

  for (uint32_t k = 0; k < m_Geometry->N_x; k++)
  {
    Index = j * m_Geometry->N_x + k;
    if(TempCol[Index]->count > 0)
    {
      for (uint32_t i = 0; i < m_Geometry->N_y; i++) //slice index
      {
        for (uint32_t q = 0; q < TempCol[Index]->count; q++)
        {
          //calculating the footprint of the voxel in the t-direction
          int16_t i_theta = int16_t(floor(static_cast<float>(TempCol[Index]->index[q] / (m_Sinogram->N_r))));
          int16_t i_r = (TempCol[Index]->index[q] % (m_Sinogram->N_r));
          VoxelLineAccessCounter = 0;
          for (uint32_t i_t = VoxelLineResponse[i]->index[0]; i_t < VoxelLineResponse[i]->index[0] + VoxelLineResponse[i]->count; i_t++) //CHANGED from <= to <
          {
            ttmp = (NuisanceParams->I_0->d[i_theta]
                * (TempCol[Index]->values[q] * VoxelLineResponse[i]->values[VoxelLineAccessCounter++] * m_Geometry->Object->getValue(j, k, i)));

            Y_Est->addToValue(ttmp, i_theta, i_r, i_t);

          }
        }
      }
    }
  }
}
