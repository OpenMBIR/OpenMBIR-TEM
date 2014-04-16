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


#include "BFUpdateYSlice.h"

//-- Boost Headers for Random Numbers
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include "MBIRLib/Common/EIMTime.h"
#include "MBIRLib/Common/EIMMath.h"

/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
#define find_min_max(low, high, V)\
{\
  low = m_Neighborhood[INDEX_3(0,0,0)];\
  high = m_Neighborhood[INDEX_3(0,0,0)];\
  for(uint8_t lcv_i = 0; lcv_i < 3;++lcv_i){\
  for(uint8_t lcv_j = 0; lcv_j < 3; ++lcv_j){\
  for(uint8_t lcv_k = 0; lcv_k < 3; ++lcv_k){\
  if(m_Neighborhood[INDEX_3(lcv_i,lcv_j,lcv_k)] < low) {low = m_Neighborhood[INDEX_3(lcv_i,lcv_j,lcv_k)];}\
  if(m_Neighborhood[INDEX_3(lcv_i,lcv_j,lcv_k)] > high) {high = m_Neighborhood[INDEX_3(lcv_i,lcv_j,lcv_k)];}\
  }\
  }\
  }\
  if(m_Theta2 !=0){\
  low = (low > (V - (m_Theta1/m_Theta2)) ? (V - (m_Theta1/m_Theta2)): low);\
  high = (high < (V - (m_Theta1/m_Theta2)) ? (V - (m_Theta1/m_Theta2)): high);\
  }\
  }


// -----------------------------------------------------------------------------
// Updates a line of voxels along y-axis
// -----------------------------------------------------------------------------
BFUpdateYSlice::BFUpdateYSlice(uint16_t yStart, uint16_t yEnd,
                           GeometryPtr geometry, int16_t outerIter, int16_t innerIter,
                           SinogramPtr  sinogram,
                           std::vector<BFAMatrixCol::Pointer> &tempCol,
                           RealVolumeType::Pointer errorSino,
                           std::vector<BFAMatrixCol::Pointer> &voxelLineResponse,
                           BFForwardModel* forwardModel,
                           UInt8Image_t::Pointer mask,
                           RealImageType::Pointer magUpdateMap, //Hold the magnitude of the reconstuction along each voxel line
                           UInt8Image_t::Pointer magUpdateMask,
                           unsigned int voxelUpdateType,
                           Real_t* averageUpdate,
                           Real_t* averageMagnitudeOfRecon,
                           unsigned int zeroSkipping,
                           BFQGGMRF::BFQGGMRF_Values *qggmrf_values,
                           VoxelUpdateList::Pointer voxelUpdateList) :
  m_YStart(yStart),
  m_YEnd(yEnd),
  m_Geometry(geometry),
  m_OuterIter(outerIter),
  m_InnerIter(innerIter),
  m_Sinogram(sinogram),
  m_TempCol(tempCol),
  m_ErrorSino(errorSino),
  m_VoxelLineResponse(voxelLineResponse),
  m_ForwardModel(forwardModel),
  m_Mask(mask),
  m_MagUpdateMap(magUpdateMap),
  m_MagUpdateMask(magUpdateMask),
 // m_VoxelUpdateType(voxelUpdateType),
  m_ZeroCount(0),
  m_CurrentVoxelValue(0.0),
  m_AverageUpdate(averageUpdate),
  m_AverageMagnitudeOfRecon(averageMagnitudeOfRecon),
  m_ZeroSkipping(zeroSkipping),
  m_QggmrfValues(qggmrf_values),
  m_VoxelUpdateList(voxelUpdateList)
{
  initVariables();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BFUpdateYSlice::~BFUpdateYSlice()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BFUpdateYSlice::initVariables()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int BFUpdateYSlice::getZeroCount()
{
  return m_ZeroCount;
}

/**
  *
  */
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
tbb::task*
#else
void
#endif
BFUpdateYSlice::execute()
{

  int32_t ArraySize = m_VoxelUpdateList->numElements();

  size_t dims[3] =
  { ArraySize, 0, 0};
  Int32ArrayType::Pointer Counter = Int32ArrayType::New(dims, "Counter");

  dims[0] = 2;
  RealArrayType::Pointer Thetas = RealArrayType::New(dims, "Thetas"); //Store
  //theta1 and theta2


  for(int32_t j = 0; j < m_VoxelUpdateList->numElements(); j++)
  {

    int32_t k_new = m_VoxelUpdateList->xIdx(j);
    int32_t j_new = m_VoxelUpdateList->zIdx(j);
    int32_t Index = j_new * m_Geometry->N_x + k_new; //This index pulls out the apprppriate index corresponding to
    //the voxel line (j_new,k_new)

    //Initialize the magnitude value to zero for appropriate pixel
    //TODO: This may be un necessary as we do this prior to calling this function ?
    m_MagUpdateMap->setValue(0, j_new, k_new);

    int shouldInitNeighborhood = 0;

    //If the Amatrix has some empty columns skip the update
    if(m_TempCol[Index]->count > 0)
    {
      ++shouldInitNeighborhood;
    }

    if(shouldInitNeighborhood > 0)
    {
      Real_t UpdatedVoxelValue = 0.0;
      int32_t errorcode = -1;
      size_t Index = j_new * m_Geometry->N_x + k_new;
      Real_t low = 0.0, high = 0.0;

      for (int32_t i = m_YStart; i < m_YEnd; i++) //slice index along Y - voxel line update
      {

        //Neighborhood of (i,j,k) should be initialized to zeros each time
        ::memset(m_Neighborhood, 0, 27*sizeof(Real_t));
        ::memset(m_BoundaryFlag, 0, 27*sizeof(uint8_t));

        //For a given (i,j,k) store its 26 point neighborhood
        for (int32_t p = -1; p <= 1; p++)
        {
          for (int32_t q = -1; q <= 1; q++)
          {
            for (int32_t r = -1; r <= 1; r++)
            {
              if(i + p >= 0 && i + p < m_Geometry->N_y)
              {
                if(j_new + q >= 0 && j_new + q < m_Geometry->N_z)
                {
                  if(k_new + r >= 0 && k_new + r < m_Geometry->N_x)
                  {
                    m_Neighborhood[INDEX_3(p + 1, q + 1, r + 1)] = m_Geometry->Object->getValue(q + j_new, r + k_new, p + i);
                    m_BoundaryFlag[INDEX_3(p + 1, q + 1, r + 1)] = 1;
                  }
                  else
                  {
                    m_BoundaryFlag[INDEX_3(p + 1, q + 1, r + 1)] = 0;
                  }
                }
              }
            }
          }
        }
        m_Neighborhood[INDEX_3(1, 1, 1)] = 0.0;
        //Compute theta1 and theta2
        m_CurrentVoxelValue = m_Geometry->Object->getValue(j_new, k_new, i); //Store the present value of the voxel
        m_Theta1 = 0.0;
        m_Theta2 = 0.0;

        //Check if every neighbor of a pixel is zero and we are not at the first iteration
        bool ZSFlag = true;
        if(m_ZeroSkipping == 1)
        {
          //Zero Skipping Algorithm
          ZSFlag = true;
          if(m_CurrentVoxelValue == 0.0 && (m_InnerIter > 0 || m_OuterIter > 0))
          {
            for (uint8_t p = 0; p <= 2; p++)
            {
              for (uint8_t q = 0; q <= 2; q++)
              {
                for (uint8_t r = 0; r <= 2; r++)
                  if(m_Neighborhood[INDEX_3(p,q,r)] > 0.0)
                  {
                    ZSFlag = false;
                    break;
                  }
              }
            }
          }
          else
          {
            ZSFlag = false; //First time dont care for zero skipping
          }
        }
        else
        {
          ZSFlag = false; //do ICD on all voxels
        }

        if(ZSFlag == false) //If the voxel is to be updated
        {
          //Forward Model parameters \theta_{1} and \theta_{2} compute
          m_ForwardModel->computeTheta(Index,m_TempCol,i,m_VoxelLineResponse,m_ErrorSino,m_Sinogram,Thetas);
          m_Theta1 = Thetas->d[0];
          m_Theta2 = Thetas->d[1];

          find_min_max(low, high, m_CurrentVoxelValue);

          if(m_Theta2 < 0){
            std::cout<<"The value of theta2 is negative"<<std::endl;
          }

          //Compute prior model parameters AND Solve the 1-D optimization problem
          errorcode = 0;
          UpdatedVoxelValue = BFQGGMRF::FunctionalSubstitution(low, high, m_CurrentVoxelValue,
                                                             m_BoundaryFlag, m_Neighborhood,
                                                             m_Theta1, m_Theta2,m_QggmrfValues);
          //Positivity constraints
          if(errorcode == 0)
          {
#ifdef POSITIVITY_CONSTRAINT
            if(UpdatedVoxelValue < 0.0)
            { //Enforcing positivity constraints
              UpdatedVoxelValue = 0.0;
            }
#endif
          }

          else {
            //Need to fill in what happens in voxel update had some numerical issues
          } //TODO Print appropriate error messages for other values of error code

          /*else //TODO: Remove this condition.
             {
               if(m_Theta1 == 0 && low == 0 && high == 0)
               {
                 UpdatedVoxelValue = 0;
               }
             }*/


          m_Geometry->Object->setValue(UpdatedVoxelValue, j_new, k_new, i);
          Real_t intermediate = m_MagUpdateMap->getValue(j_new, k_new) + fabs(UpdatedVoxelValue - m_CurrentVoxelValue);
          m_MagUpdateMap->setValue(intermediate, j_new, k_new);

#if ROI
          if(m_Mask->getValue(j_new, k_new) == 1)
          { //Stopping criteria variables for "Full update ICD" algorithm
            *m_AverageUpdate += fabs(UpdatedVoxelValue - m_CurrentVoxelValue);
            *m_AverageMagnitudeOfRecon += fabs(m_CurrentVoxelValue); //computing the percentage update =(Change in mag/Initial magnitude)
          }
#endif //ROI
          //Update the ErrorSinogram and Bragg selector
          m_ForwardModel->updateErrorSinogram(UpdatedVoxelValue - m_CurrentVoxelValue, Index, m_TempCol, i, m_VoxelLineResponse, m_ErrorSino, m_Sinogram);
        }
        else
        {
          m_ZeroCount++;
        }

      }
    }
    else
    {
      continue;
    }


  }

#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
  return NULL;
#endif
}

