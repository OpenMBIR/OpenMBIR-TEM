/*
 * ForwardProject.h
 *
 *  Created on: Dec 9, 2011
 *      Author: mjackson
 */

#ifndef FORWARDPROJECT_H_
#define FORWARDPROJECT_H_

#include <iostream>


#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class ForwardProject
{
  public:
    ForwardProject(Sinogram* sinogram, Geometry* geometry,
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
    {}
    virtual ~ForwardProject(){};
    void operator()() const
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
  private:
    Sinogram* m_Sinogram;
    Geometry* m_Geometry;
    AMatrixCol*** TempCol;
    AMatrixCol* VoxelLineResponse;
    RealVolumeType::Pointer Y_Est;
    ScaleOffsetParams* NuisanceParams;
    uint16_t m_Tilt;
};

#endif /* FORWARDPROJECT_H_ */
