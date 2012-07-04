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


#include "ComputeInitialOffsets.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeInitialOffsets::ComputeInitialOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeInitialOffsets::~ComputeInitialOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//Initialize offsets via least squares fit
void ComputeInitialOffsets::execute()
{
  // If an error occurs, clean up any memory, call "setErrorCondition(-1)" and
  // also setErrorMessage("Something went wrong"); and then return

  notify("GainsOffsetsCalculation Starting", 0, UpdateProgressMessage);
  SinogramPtr sinogram = getSinogram(); //This I assume some how gets the sinogram as it stands now
  TomoInputsPtr inputs = getTomoInputs(); //This gets the input files

  //The normalization and offset parameters for the views

//  size_t dims[3] =
//  { sinogram->N_theta, 0, 0 };
//  sinogram->InitialGain = RealArrayType::New(dims);
//  sinogram->InitialGain->setName("sinogram->InitialGain");
//  sinogram->InitialOffset = RealArrayType::New(dims);
//  sinogram->InitialOffset->setName("sinogram->InitialOffset");
//  sinogram->InitialVariance = RealArrayType::New(dims);
//  sinogram->InitialVariance->setName("sinogram->InitialVariance");
  //Form the average Gain per view
  std::vector<Real_t> AverageGain(sinogram->N_theta, 0);
  std::vector<Real_t> TargetGain(sinogram->N_theta, 0);
//	std::vector<std::vector<DATA_TYPE>>LS_Matrix(sinogram->N_theta, std::vector<DATA_TYPE> (2));
//  Real_t** LS_Matrix = (Real_t**)get_img(2, sinogram->N_theta, sizeof(Real_t));

  size_t dims[2] = {sinogram->N_theta, 2};
  RealImageType::Pointer LS_Matrix = RealImageType::New(dims, "LS_Matrix");


  Real_t LS_Estimates[2] = {0, 0};

  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
  {
    Real_t sum = 0;
    for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
      {
        sum += sinogram->counts->getValue(i_theta, i_r, i_t);
      }
    }
    sum /= (sinogram->N_r * sinogram->N_t);
    AverageGain[i_theta] = sum;
    TargetGain[i_theta] = 1.0 / cos(sinogram->angles[i_theta] * M_PI / 180); //Set to 1/cos(tilt_angle)
    LS_Matrix->setValue(TargetGain[i_theta], i_theta, 0);
    LS_Matrix->setValue(1, i_theta, 1);
   // LS_Matrix[i_theta][1] = 1;
  }

#if 0
  std::cout << "Average Gains" << std::endl;
  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
  {
    std::cout << AverageGain[i_theta] << std::endl;
  }

  std::cout << "Target Gains" << std::endl;
  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
  {
    std::cout << TargetGain[i_theta] << std::endl;
  }
#endif

//	DATA_TYPE max_elt=*max_element(AverageGain.begin(), AverageGain.end());
//	DATA_TYPE min_elt=*min_element(AverageGain.begin(), AverageGain.end());

  //Compute A^t * A
  Real_t ProdMatrix[2][2];

  for (uint8_t i = 0; i < 2; i++)
  {
    for (uint8_t j = 0; j < 2; j++)
    {
      ProdMatrix[i][j] = 0;
      for (uint8_t k = 0; k < sinogram->N_theta; k++)
      {
        ProdMatrix[i][j] += (LS_Matrix->getValue(k, i) * LS_Matrix->getValue(k, j));
      }
    }
  }

  //Compute T1=(A^t*A)^-1
  Real_t determinant = ProdMatrix[0][0] * ProdMatrix[1][1] - (ProdMatrix[0][1] * ProdMatrix[0][1]);
  Real_t InverseProdMatrix[2][2];

  InverseProdMatrix[0][0] = ProdMatrix[1][1] / determinant;
  InverseProdMatrix[1][1] = ProdMatrix[0][0] / determinant;
  InverseProdMatrix[0][1] = -ProdMatrix[1][0] / determinant;
  InverseProdMatrix[1][0] = -ProdMatrix[0][1] / determinant;

  //Compute T2=A^t*Y
  std::vector<Real_t> TempVector(2, 0);
  for (uint8_t j = 0; j < 2; j++)
  {
    for (uint8_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      TempVector[j] += LS_Matrix->getValue(i_theta, j) * AverageGain[i_theta];
    }
  }


  //Compute T1*T2 (2 X 2 matrix with 2 X 1 vector)
  for (uint8_t i = 0; i < 2; i++)
  {
    for (uint8_t k = 0; k < 2; k++)
    {
      LS_Estimates[i] += InverseProdMatrix[i][k] * TempVector[k];
    }
  }

  std::stringstream ss;
  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
  {
    sinogram->InitialOffset->d[i_theta] = LS_Estimates[1];
    ss << "Tilt: " << i_theta << "  Offset: " << sinogram->InitialOffset->d[i_theta] << std::endl;
  }
  if (getVeryVerbose()) {
    std::cout << ss.str() << std::endl;
  }
  setErrorCondition(0);
  setErrorMessage("");
  notify("Done ComputeInitialOffsets", 0, UpdateProgressMessage);
}
