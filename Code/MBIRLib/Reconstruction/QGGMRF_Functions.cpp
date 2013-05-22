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



#include "QGGMRF_Functions.h"

#include "MBIRLib/Common/EIMMath.h"



namespace Detail {
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  inline double Clip(double x, double a, double b)
  {
    return (x < a) ? a : ((x > b) ? b:x);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  inline int16_t mod(int16_t a,int16_t b)
  {
    int16_t temp;
    temp=a%b;
    if(temp < 0)
      return temp + b;
    else {
      return temp;
    }

  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  inline double Minimum(double a, double b)
  {
    return (a < b ? a: b);
  }
}

namespace QGGMRF {

void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF_Values* qggmrf_values)
{
    qggmrf_values->MRF_P = 2;
    qggmrf_values->MRF_Q = tomoInputs->p;
    qggmrf_values->MRF_C = 0.001;
    qggmrf_values->MRF_ALPHA = 1.5;
    qggmrf_values->SIGMA_X_P = pow(tomoInputs->SigmaX, qggmrf_values->MRF_P);
    qggmrf_values->SIGMA_X_P_Q = pow(tomoInputs->SigmaX, (qggmrf_values->MRF_P - qggmrf_values->MRF_Q));
    qggmrf_values->SIGMA_X_Q = pow(tomoInputs->SigmaX, qggmrf_values->MRF_Q);
	qggmrf_values->gamma = 5;
}
	
	
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t updatePriorModel(Real_t NewSigmaX,QGGMRF_Values* qggmrf_values)
{
    Real_t SigmaX = qggmrf_values->gamma*pow(NewSigmaX*qggmrf_values->MRF_Q,1.0/qggmrf_values->MRF_Q);
	qggmrf_values->SIGMA_X_P = pow(SigmaX, qggmrf_values->MRF_P);
	qggmrf_values->SIGMA_X_P_Q = pow(SigmaX, (qggmrf_values->MRF_P - qggmrf_values->MRF_Q));
	qggmrf_values->SIGMA_X_Q = pow(SigmaX, qggmrf_values->MRF_Q);
	return SigmaX;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t Value(Real_t delta, QGGMRF::QGGMRF_Values* qggmrf_values)
{
  return ((pow(fabs(delta), qggmrf_values->MRF_P) / qggmrf_values->SIGMA_X_P) / (qggmrf_values->MRF_C + pow(fabs(delta), qggmrf_values->MRF_P - qggmrf_values->MRF_Q) / qggmrf_values->SIGMA_X_P_Q));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t Derivative(Real_t delta, QGGMRF::QGGMRF_Values* qggmrf_values)
{
  Real_t temp1, temp2, temp3;
  temp1 = pow(fabs(delta), qggmrf_values->MRF_P - qggmrf_values->MRF_Q) / (qggmrf_values->SIGMA_X_P_Q);
  temp2 = pow(fabs(delta), qggmrf_values->MRF_P - 1);
  temp3 = qggmrf_values->MRF_C + temp1;
  if(delta < 0)
  {
    return ((-1 * temp2 / (temp3 * qggmrf_values->SIGMA_X_P)) * (qggmrf_values->MRF_P - ((qggmrf_values->MRF_P - qggmrf_values->MRF_Q) * temp1) / (temp3)));
  }
  else
  {
    return ((temp2 / (temp3 * qggmrf_values->SIGMA_X_P)) * (qggmrf_values->MRF_P - ((qggmrf_values->MRF_P - qggmrf_values->MRF_Q) * temp1) / (temp3)));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t SecondDerivative(Real_t delta, QGGMRF::QGGMRF_Values* qggmrf_values)
{
  return qggmrf_values->MRF_P / (qggmrf_values->SIGMA_X_P * qggmrf_values->MRF_C);
}

// -----------------------------------------------------------------------------
//       Real_t QGGMRF_Params[26][3];
// Function to compute parameters of thesurrogate function
// -----------------------------------------------------------------------------
void ComputeParameters(Real_t umin, Real_t umax, Real_t RefValue,
                       uint8_t* BOUNDARYFLAG, Real_t* NEIGHBORHOOD,
                       QGGMRF::QGGMRF_Values* qggmrf_values,
                       Real_t* QGGMRF_Params)
{
  Real_t Delta0;
  uint8_t i, j, k, count = 0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        if((i != 1 || j != 1 || k != 1) && BOUNDARYFLAG[INDEX_3(i,j,k)] == 1)
        {
          Delta0 = RefValue - NEIGHBORHOOD[INDEX_3(i,j,k)];
          if(Delta0 != 0)
          {
            QGGMRF_Params[count*3 + 0] = QGGMRF::Derivative(Delta0, qggmrf_values) / (Delta0);
          }
          else
          {
            QGGMRF_Params[count*3 + 0] = QGGMRF::SecondDerivative(0, qggmrf_values);
          }
          count++;
        }
      }
    }
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t FunctionalSubstitution(Real_t umin, Real_t umax, Real_t currentVoxelValue,
                              uint8_t* BOUNDARYFLAG, Real_t* FILTER, Real_t* NEIGHBORHOOD,
                              Real_t THETA1, Real_t THETA2,
                              QGGMRF::QGGMRF_Values* qggmrf_values)
{
  Real_t u, temp1 = 0, temp2 = 0, temp_const, RefValue = 0;
  uint8_t i, j, k, count = 0;
#ifdef POSITIVITY_CONSTRAINT
  if(umin < 0)
  umin =0;
#endif //Positivity
  RefValue = currentVoxelValue;
  //Need to Loop this for multiple iterations of substitute function
  Real_t QGGMRF_Params[26*3];
  for (uint8_t qggmrf_iter = 0; qggmrf_iter < QGGMRF::QGGMRF_ITER; qggmrf_iter++)
  {
    QGGMRF::ComputeParameters(umin, umax, RefValue, BOUNDARYFLAG, NEIGHBORHOOD,
                              qggmrf_values, QGGMRF_Params);
	  count = 0;// cycles through the list of parameters of the quadratic surrogate
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        for (k = 0; k < 3; k++)
        {
          if((i != 1 || j != 1 || k != 1) && BOUNDARYFLAG[INDEX_3(i,j,k)] == 1)
          {
            temp_const = FILTER[INDEX_3(i,j,k)] * QGGMRF_Params[count*3 + 0];
            temp1 += temp_const * NEIGHBORHOOD[INDEX_3(i,j,k)];
            temp2 += temp_const;
            count++;
          }
        }
      }
    }
    u = (temp1 + (THETA2 * currentVoxelValue) - THETA1) / (temp2 + THETA2);

    if(qggmrf_iter < QGGMRF::QGGMRF_ITER - 1)
    {
      RefValue = Detail::Clip(RefValue + qggmrf_values->MRF_ALPHA * (u - RefValue), umin, umax);
    }
    else
    {
      RefValue = Detail::Clip(RefValue + qggmrf_values->MRF_ALPHA * (u - RefValue), umin, umax);
    }
  }

  return RefValue;
}


} /* End Namespace */


