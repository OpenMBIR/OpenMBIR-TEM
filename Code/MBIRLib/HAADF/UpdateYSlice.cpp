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


#include "UpdateYSlice.h"

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
//
// -----------------------------------------------------------------------------
UpdateYSlice::UpdateYSlice(uint16_t yStart, uint16_t yEnd,
                           GeometryPtr geometry, int16_t outerIter, int16_t innerIter,
                           SinogramPtr  sinogram,
                           std::vector<HAADFAMatrixCol::Pointer> &tempCol,
                           RealVolumeType::Pointer errorSino,
                           std::vector<HAADFAMatrixCol::Pointer> &voxelLineResponse,
                           HAADFForwardModel* forwardModel,
                           UInt8Image_t::Pointer mask,
                           RealImageType::Pointer magUpdateMap, //Hold the magnitude of the reconstuction along each voxel line
                           UInt8Image_t::Pointer magUpdateMask,
                           unsigned int voxelUpdateType,
                           Real_t nh_Threshold,
                           Real_t* averageUpdate,
                           Real_t* averageMagnitudeOfRecon,
                           unsigned int zeroSkipping,
                           QGGMRF::QGGMRF_Values *qggmrf_values) :
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
m_VoxelUpdateType(voxelUpdateType),
m_NH_Threshold(nh_Threshold),
m_ZeroCount(0),
m_CurrentVoxelValue(0.0),
m_AverageUpdate(averageUpdate),
m_AverageMagnitudeOfRecon(averageMagnitudeOfRecon),
m_ZeroSkipping(zeroSkipping),
m_QggmrfValues(qggmrf_values)
{
  initVariables();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
UpdateYSlice::~UpdateYSlice()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void UpdateYSlice::initVariables()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int UpdateYSlice::getZeroCount()
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
   UpdateYSlice::execute()
   {

     const uint32_t rangeMin = 0;
     const uint32_t rangeMax = std::numeric_limits<uint32_t>::max();
     typedef boost::uniform_int<uint32_t> NumberDistribution;
     typedef boost::mt19937 RandomNumberGenerator;
     typedef boost::variate_generator<RandomNumberGenerator&,
                                      NumberDistribution> Generator;
     NumberDistribution distribution(rangeMin, rangeMax);
     RandomNumberGenerator generator;
     Generator numberGenerator(generator, distribution);
     boost::uint32_t arg = static_cast<boost::uint32_t>(EIMTOMO_getMilliSeconds());
     generator.seed(arg); // seed with the current time
	   
     int32_t ArraySize = m_Geometry->N_x * m_Geometry->N_z;
     size_t dims[3] =
     { ArraySize, 0, 0};
     Int32ArrayType::Pointer Counter = Int32ArrayType::New(dims, "Counter");

	dims[0] = 2;
	RealArrayType::Pointer Thetas = RealArrayType::New(dims, "Thetas");
	   
     dims[0] = m_Geometry->N_z;
     dims[1] = m_Geometry->N_x;
     dims[2] = 0;
     UInt8Image_t::Pointer m_VisitCount = UInt8Image_t::New(dims, "VisitCount");

     for (int32_t j_new = 0; j_new < ArraySize; j_new++)
     {
       Counter->d[j_new] = j_new;
     }
     uint32_t NumVoxelsToUpdate = 0;
     for (int32_t j = 0; j < m_Geometry->N_z; j++)
     {
       for (int32_t k = 0; k < m_Geometry->N_x; k++)
       {
         if(m_VoxelUpdateType == VoxelUpdateType::NonHomogeniousUpdate)
         {
           if(m_MagUpdateMap->getValue(j, k) > m_NH_Threshold)
           {
             m_MagUpdateMask->setValue(1, j, k);
             m_MagUpdateMap->setValue(0, j, k);
             NumVoxelsToUpdate++;
           }
           else
           {
             m_MagUpdateMask->setValue(0, j, k);
           }
         }
         else if(m_VoxelUpdateType == VoxelUpdateType::HomogeniousUpdate)
         {
           m_MagUpdateMap->setValue(0, j, k);
           NumVoxelsToUpdate++;
         }
         else if(m_VoxelUpdateType == VoxelUpdateType::RegularRandomOrderUpdate)
         {
           m_MagUpdateMap->setValue(0, j, k);
           NumVoxelsToUpdate++;
         }
         m_VisitCount->setValue(0, j, k);
       }
     }
     //   std::cout << "    " << "Number of voxel lines to update: " << NumVoxelsToUpdate << std::endl;
     for (int32_t j = 0; j < m_Geometry->N_z; j++) //Row index
     {
       for (int32_t k = 0; k < m_Geometry->N_x; k++) //Column index
       {
         uint32_t Index = numberGenerator() % ArraySize;
		 int32_t k_new = Counter->d[Index] % m_Geometry->N_x;
		 int32_t j_new = Counter->d[Index] / m_Geometry->N_x;
         Counter->d[Index] = Counter->d[ArraySize - 1];
         m_VisitCount->setValue(1, j_new, k_new);
         ArraySize--;
         Index = j_new * m_Geometry->N_x + k_new; //This index pulls out the apprppriate index corresponding to
         //the voxel line (j_new,k_new)
         int shouldInitNeighborhood = 0;

         if(m_VoxelUpdateType == VoxelUpdateType::NonHomogeniousUpdate
             && m_MagUpdateMask->getValue(j_new, k_new) == 1
             && m_TempCol[Index]->count > 0)
         {
           ++shouldInitNeighborhood;
         }
         if(m_VoxelUpdateType == VoxelUpdateType::HomogeniousUpdate
             && m_TempCol[Index]->count > 0)
         {
           ++shouldInitNeighborhood;
         }
         if(m_VoxelUpdateType == VoxelUpdateType::RegularRandomOrderUpdate
             && m_TempCol[Index]->count > 0)
         {
           ++shouldInitNeighborhood;
         }

         if(shouldInitNeighborhood > 0)
         //After this should ideally call UpdateVoxelLine(j_new,k_new) ie put everything in this "if" inside a method called UpdateVoxelLine
         {
           //   std::cout << "UpdateYSlice- YStart: " << m_YStart << "  YEnd: " << m_YEnd << std::endl;
           Real_t UpdatedVoxelValue = 0.0;
           int32_t errorcode = -1;
           size_t Index = j_new * m_Geometry->N_x + k_new;
           Real_t low = 0.0, high = 0.0;

           for (int32_t i = m_YStart; i < m_YEnd; i++) //slice index
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

             if(ZSFlag == false)
             {
			  //Forward Model parameters \theta_{1} and \theta_{2} compute
			  m_ForwardModel->computeTheta(Index,m_TempCol,i,m_VoxelLineResponse,m_ErrorSino,m_Sinogram,Thetas);
			  m_Theta1 = Thetas->d[0];
			  m_Theta2 = Thetas->d[1];
			  find_min_max(low, high, m_CurrentVoxelValue);
			
				 if(m_Theta2 < 0){
				//std::cout<<"The value of theta2 is negative"<<std::endl;
				 }
			  	 
               //Compute prior model parameters AND Solve the 1-D optimization problem
               errorcode = 0;
			   UpdatedVoxelValue = QGGMRF::FunctionalSubstitution(low, high, m_CurrentVoxelValue,
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
               else //TODO Remove this condition.
               {
                 if(m_Theta1 == 0 && low == 0 && high == 0)
                 {
                   UpdatedVoxelValue = 0;
                 }
               }

               //TODO Print appropriate error messages for other values of error code
               m_Geometry->Object->setValue(UpdatedVoxelValue, j_new, k_new, i);
               Real_t intermediate = m_MagUpdateMap->getValue(j_new, k_new) + fabs(UpdatedVoxelValue - m_CurrentVoxelValue);
               m_MagUpdateMap->setValue(intermediate, j_new, k_new);

#if ROI
               //if(Mask->d[j_new][k_new] == 1)
               if(m_Mask->getValue(j_new, k_new) == 1)
               {
                 *m_AverageUpdate += fabs(UpdatedVoxelValue - m_CurrentVoxelValue);
                 *m_AverageMagnitudeOfRecon += fabs(m_CurrentVoxelValue); //computing the percentage update =(Change in mag/Initial magnitude)
               }
#endif //ROI
               
               //Update the ErrorSinogram
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
     }

     for (int j = 0; j < m_Geometry->N_z; j++)
     { //Row index
       for (int k = 0; k < m_Geometry->N_x; k++)
       { //Column index
         if(m_VisitCount->getValue(j, k) == 0)
         {
           printf("Pixel (%d %d) not visited\n", j, k);
         }
       }
     }
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
     return NULL;
#endif
   }

