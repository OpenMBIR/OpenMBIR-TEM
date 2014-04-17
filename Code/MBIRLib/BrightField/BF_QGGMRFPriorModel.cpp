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
#include "BF_QGGMRFPriorModel.h"

#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/BrightField/BFConstants.h"


namespace QGGMRF {

  // -----------------------------------------------------------------------------
  //Contains all the prior model related functions
  // -----------------------------------------------------------------------------
  void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF::QGGMRF_Values* qggmrf_values, Real_t* filter)
  {
#ifdef EIMTOMO_USE_QGGMRF
    qggmrf_values->MRF_P = 2;
    qggmrf_values->MRF_Q = tomoInputs->p;
    qggmrf_values->MRF_C = 0.01;
    qggmrf_values->MRF_ALPHA = 1.5;
    qggmrf_values->SIGMA_X_P = pow(tomoInputs->SigmaX, qggmrf_values->MRF_P);
    qggmrf_values->SIGMA_X_P_Q = pow(tomoInputs->SigmaX, (qggmrf_values->MRF_P - qggmrf_values->MRF_Q));
    qggmrf_values->SIGMA_X_Q = pow(tomoInputs->SigmaX, qggmrf_values->MRF_Q);
#else
    MRF_P = tomoInputs->p;
    SIGMA_X_P = pow(tomoInputs->SigmaX,MRF_P);
#endif //QGGMRF
    if(NULL != filter)
    {
      initializeFilter(filter);
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void initializeFilter(Real_t* filter)
  {
    filter[INDEX_3(0,0,0)] = 0.0302;
    filter[INDEX_3(0,0,1)] = 0.0370;
    filter[INDEX_3(0,0,2)] = 0.0302;
    filter[INDEX_3(0,1,0)] = 0.0370;
    filter[INDEX_3(0,1,1)] = 0.0523;
    filter[INDEX_3(0,1,2)] = 0.0370;
    filter[INDEX_3(0,2,0)] = 0.0302;
    filter[INDEX_3(0,2,1)] = 0.0370;
    filter[INDEX_3(0,2,2)] = 0.0302;

    filter[INDEX_3(1,0,0)] = 0.0370;
    filter[INDEX_3(1,0,1)] = 0.0523;
    filter[INDEX_3(1,0,2)] = 0.0370;
    filter[INDEX_3(1,1,0)] = 0.0523;
    filter[INDEX_3(1,1,1)] = 0.0000;
    filter[INDEX_3(1,1,2)] = 0.0523;
    filter[INDEX_3(1,2,0)] = 0.0370;
    filter[INDEX_3(1,2,1)] = 0.0523;
    filter[INDEX_3(1,2,2)] = 0.0370;

    filter[INDEX_3(2,0,0)] = 0.0302;
    filter[INDEX_3(2,0,1)] = 0.0370;
    filter[INDEX_3(2,0,2)] = 0.0302;
    filter[INDEX_3(2,1,0)] = 0.0370;
    filter[INDEX_3(2,1,1)] = 0.0523;
    filter[INDEX_3(2,1,2)] = 0.0370;
    filter[INDEX_3(2,2,0)] = 0.0302;
    filter[INDEX_3(2,2,1)] = 0.0370;
    filter[INDEX_3(2,2,2)] = 0.0302;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  Real_t PriorModelCost(GeometryPtr geometry, QGGMRF_Values* qggmrf_values, Real_t* filter)
  {
    Real_t cost = 0;
    Real_t temp = 0;
    Real_t delta=0;

#ifndef EIMTOMO_USE_QGGMRF
    for (int16_t i = 0; i < geometry->N_z; i++)
      for (int16_t j = 0; j < geometry->N_x; j++)
        for (int16_t k = 0; k < geometry->N_y; k++)
        {

          if(k + 1 < geometry->N_y) temp += k_Filter[2][1][1] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i][j][k + 1]), MRF_P);

          if(j + 1 < geometry->N_x)
          {
            if(k - 1 >= 0) temp += k_Filter[0][1][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i][j + 1][k - 1]), MRF_P);

            temp += k_Filter[1][1][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i][j + 1][k]), MRF_P);

            if(k + 1 < geometry->N_y) temp += k_Filter[2][1][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i][j + 1][k + 1]), MRF_P);

          }

          if(i + 1 < geometry->N_z)
          {

            if(j - 1 >= 0) temp += k_Filter[1][2][0] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j - 1][k]), MRF_P);

            temp += k_Filter[1][2][1] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j][k]), MRF_P);

            if(j + 1 < geometry->N_x) temp += k_Filter[1][2][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j + 1][k]), MRF_P);

            if(j - 1 >= 0)
            {
              if(k - 1 >= 0) temp += k_Filter[0][2][0] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j - 1][k - 1]), MRF_P);

              if(k + 1 < geometry->N_y) temp += k_Filter[2][2][0] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j - 1][k + 1]), MRF_P);

            }

            if(k - 1 >= 0) temp += k_Filter[0][2][1] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j][k - 1]), MRF_P);

            if(j + 1 < geometry->N_x)
            {
              if(k - 1 >= 0) temp += k_Filter[0][2][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j + 1][k - 1]), MRF_P);

              if(k + 1 < geometry->N_y) temp += k_Filter[2][2][2] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j + 1][k + 1]), MRF_P);
            }

            if(k + 1 < geometry->N_y) temp += k_Filter[2][2][1] * pow(fabs(geometry->Object->d[i][j][k] - geometry->Object->d[i + 1][j][k + 1]), MRF_P);
          }
        }
    cost += (temp / (MRF_P * SIGMA_X_P));
#else

    for (int16_t i = 0; i < geometry->N_z; i++)
    {
      for (int16_t j = 0; j < geometry->N_x; j++)
      {
        for (int16_t k = 0; k < geometry->N_y; k++)
        {

          if(k + 1 < geometry->N_y)
          {
            delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i, j, k + 1);
            temp += filter[INDEX_3(2,1,1)] * QGGMRF::Value(delta, qggmrf_values);

          }

          if(j + 1 < geometry->N_x)
          {
            if(k - 1 >= 0)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i, j + 1, k - 1);
              temp += filter[INDEX_3(0,1,2)] * QGGMRF::Value(delta, qggmrf_values);
            }

            delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i, j + 1, k);
            temp += filter[INDEX_3(1,1,2)] * QGGMRF::Value(delta, qggmrf_values);

            if(k + 1 < geometry->N_y)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i, j + 1, k + 1);
              temp += filter[INDEX_3(2,1,2)] * QGGMRF::Value(delta, qggmrf_values);
            }

          }

          if(i + 1 < geometry->N_z)
          {

            if(j - 1 >= 0)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j - 1, k);
              temp += filter[INDEX_3(1,2,0)] * QGGMRF::Value(delta, qggmrf_values);
            }

            delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j, k);
            temp += filter[INDEX_3(1,2,1)] * QGGMRF::Value(delta, qggmrf_values);

            if(j + 1 < geometry->N_x)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j + 1, k);
              temp += filter[INDEX_3(1,2,2)] * QGGMRF::Value(delta, qggmrf_values);
            }

            if(j - 1 >= 0)
            {
              if(k - 1 >= 0)
              {
                delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j - 1, k - 1);
                temp += filter[INDEX_3(0,2,0)] * QGGMRF::Value(delta, qggmrf_values);
              }

              if(k + 1 < geometry->N_y)
              {
                delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j - 1, k + 1);
                temp += filter[INDEX_3(2,2,0)] * QGGMRF::Value(delta, qggmrf_values);
              }

            }

            if(k - 1 >= 0)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j, k - 1);
              temp += filter[INDEX_3(0,2,1)] * QGGMRF::Value(delta, qggmrf_values);
            }

            if(j + 1 < geometry->N_x)
            {
              if(k - 1 >= 0)
              {
                delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j + 1, k - 1);
                temp += filter[INDEX_3(0,2,2)] * QGGMRF::Value(delta, qggmrf_values);
              }

              if(k + 1 < geometry->N_y)
              {
                delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j + 1, k + 1);
                temp += filter[INDEX_3(2,2,2)] * QGGMRF::Value(delta, qggmrf_values);
              }
            }

            if(k + 1 < geometry->N_y)
            {
              delta = geometry->Object->getValue(i, j, k) - geometry->Object->getValue(i + 1, j, k + 1);
              temp += filter[INDEX_3(2,2,1)] * QGGMRF::Value(delta, qggmrf_values);
            }
          }
        }
      }
    }
    cost += (temp);
#endif //QGGMRF

    return cost;
  }


} /* End Namespace */


