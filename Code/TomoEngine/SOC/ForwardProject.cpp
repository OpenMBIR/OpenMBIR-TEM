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

#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ForwardProject::ForwardProject(Sinogram* sinogram, Geometry* geometry,
                               AMatrixCol*** tempCol, AMatrixCol* voxelLineResponse,
                               RealVolumeType::Pointer yEst,
                               ScaleOffsetParams* nuisanceParams,
                               uint16_t tilt) :
                  m_Sinogram(sinogram),
                  m_Geometry(geometry),
                  TempCol(tempCol),
                  VoxelLineResponse(voxelLineResponse),
                  Y_Est(yEst),
                  NuisanceParams(nuisanceParams),
                  m_Tilt(tilt)
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
  std::cout << "Forward projecting Tilt " << m_Tilt << "/" << m_Geometry->N_z << std::endl;

  uint16_t j = m_Tilt;
  uint16_t VoxelLineAccessCounter;
  for (uint16_t k = 0; k < m_Geometry->N_x; k++)
  {
    if(TempCol[j][k]->count > 0)
    {
      for (uint16_t i = 0; i < m_Geometry->N_y; i++) //slice index
      {
        for (uint32_t q = 0; q < TempCol[j][k]->count; q++)
        {
          //calculating the footprint of the voxel in the t-direction
          int16_t i_theta = int16_t(floor(static_cast<float>(TempCol[j][k]->index[q] / (m_Sinogram->N_r))));
          int16_t i_r = (TempCol[j][k]->index[q] % (m_Sinogram->N_r));

          VoxelLineAccessCounter = 0;
          for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++) //CHANGED from <= to <
          {
            Y_Est->d[i_theta][i_r][i_t] += (NuisanceParams->I_0->d[i_theta] * (TempCol[j][k]->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter++] * m_Geometry->Object->d[j][k][i]));

          }
        }
      }
    }
  }
}
