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

#include "VoxelUpdate.h"

#include <stdio.h>

#include "TomoEngine/Common/EIMTime.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/DerivOfCostFunc.hpp"



#define START startm = EIMTOMO_getMilliSeconds();
#define STOP stopm = EIMTOMO_getMilliSeconds();
#define PRINTTIME printf( "%6.3f seconds used by the processor.\n", ((double)stopm-startm)/1000.0);



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdate::VoxelUpdate() :
    TomoFilter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdate::~VoxelUpdate()
{
}

/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void VoxelUpdate::minMax(DATA_TYPE *low,DATA_TYPE *high)
{
  DATA_TYPE THETA1 = m_VoxelUpdateValues->THETA1;
  DATA_TYPE THETA2 = m_VoxelUpdateValues->THETA2;
  DATA_TYPE V = m_VoxelUpdateValues->V;

  *low=m_VoxelUpdateValues->NEIGHBORHOOD[0][0][0];
  *high=m_VoxelUpdateValues->NEIGHBORHOOD[0][0][0];

  for(int32_t i = 0; i < 3;i++)
  {
    for(int32_t j=0; j < 3; j++)
    {
      for(int32_t k = 0; k < 3; k++)
      {
        //  if(NEIGHBORHOOD[i][j][k] != 0)
        //  printf("%lf ", NEIGHBORHOOD[i][j][k]);

        if(m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k] < *low) {
          *low = m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k];}
        if(m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k] > *high){
          *high=m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k];}
      }
      //  printf("\n");
    }
  }


  if(THETA2 !=0)
  {
  *low = (*low > (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *low);
  *high = (*high < (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *high);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdate::execute()
{
  // Copy from the structure into local variables purely for short syntax
  int16_t Iter = m_VoxelUpdateValues->Iter;
  uint16_t* cost_counter = m_VoxelUpdateValues->cost_counter;
  Int32ArrayType::Pointer Counter = m_VoxelUpdateValues->Counter;
  UInt8ImageType::Pointer VisitCount = m_VoxelUpdateValues->VisitCount;
  AMatrixCol*** TempCol = m_VoxelUpdateValues->TempCol;
  AMatrixCol* VoxelLineResponse = m_VoxelUpdateValues->VoxelLineResponse;
  ScaleOffsetParams* NuisanceParams = m_VoxelUpdateValues->NuisanceParams;
  RealVolumeType::Pointer ErrorSino = m_VoxelUpdateValues->ErrorSino;
  RealVolumeType::Pointer Weight = m_VoxelUpdateValues->Weight;
  UInt8ImageType::Pointer Mask = m_VoxelUpdateValues->Mask;
  RealArrayType::Pointer cost = m_VoxelUpdateValues->cost;
  RNGVars* RandomNumber = m_VoxelUpdateValues->RandomNumber;
 // uint8_t BOUNDARYFLAG[3][3][3] = m_VoxelUpdateValues->BOUNDARYFLAG;
 // DATA_TYPE FILTER[3][3][3] = m_VoxelUpdateValues->FILTER;
  DATA_TYPE THETA1 = m_VoxelUpdateValues->THETA1;
  DATA_TYPE THETA2 = m_VoxelUpdateValues->THETA2;
 // DATA_TYPE NEIGHBORHOOD[3][3][3] = m_VoxelUpdateValues->NEIGHBORHOOD;
  DATA_TYPE V = m_VoxelUpdateValues->V;
  DATA_TYPE MRF_P = m_VoxelUpdateValues->MRF_P;
  DATA_TYPE SIGMA_X_P = m_VoxelUpdateValues->SIGMA_X_P;
  FILE* Fp2 = m_VoxelUpdateValues->Fp2;



  // These are local variables
  int32_t ArraySize;
  int32_t Index;
  int32_t j_new, k_new;
  uint16_t VoxelLineAccessCounter;
  DATA_TYPE low, high;
  DATA_TYPE accuracy = 1e-7; //This is the rooting accuracy for x
  int32_t errorcode = -1;
  AMatrixCol* TempMemBlock;
  uint64_t startm = 0;
  uint64_t stopm = 0;

  // Get pointers to the Sinogram and Geometry
  Sinogram* sinogram = getSinogram();
  Geometry* geometry = getGeometry();

#ifdef ROI
  //variables used to stop the process
  DATA_TYPE AverageUpdate;
  DATA_TYPE AverageMagnitudeOfRecon;
#endif

  high = (DATA_TYPE)INT64_MAX;
  low = (DATA_TYPE)INT64_MIN;

  //printf("Iter %d\n",Iter);
#ifdef ROI
  AverageUpdate = 0;
  AverageMagnitudeOfRecon = 0;
#endif

#ifdef RANDOM_ORDER_UPDATES
  ArraySize = geometry->N_x * geometry->N_z;
  for (int32_t j = 0; j < ArraySize; j++)
  {
    Counter->d[j] = j;
  }

  for (uint16_t j = 0; j < geometry->N_z; j++)
  {
    for (uint16_t k = 0; k < geometry->N_x; k++)
    {
      VisitCount->d[j][k] = 0;
    }
  }
#endif

  START;
  for (uint16_t j = 0; j < geometry->N_z; j++) //Row index
  {
    for (uint16_t k = 0; k < geometry->N_x; k++) //Column index
    {
#ifdef RANDOM_ORDER_UPDATES
      //RandomNumber=init_genrand(Iter);
      Index = (genrand_int31(RandomNumber)) % ArraySize;
      k_new = Counter->d[Index] % geometry->N_x;
      j_new = Counter->d[Index] / geometry->N_x;
      //memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index-1));
      //TODO: Instead just swap the value in Index with the one in ArraySize
      Counter->d[Index] = Counter->d[ArraySize - 1];
      VisitCount->d[j_new][k_new] = 1;
      ArraySize--;
#else
      j_new = j;
      k_new = k;
#endif //Random order updates
      TempMemBlock = TempCol[j_new][k_new]; //Remove this
      if(TempMemBlock->count > 0)
      {
        int32_t Idx = 0;
        for (int32_t i = 0; i < geometry->N_y; i++) //slice index
        {
          //Neighborhood of (i,j,k) should be initialized to zeros each time
          for (int32_t p = 0; p <= 2; p++)
          {
            for (int32_t q = 0; q <= 2; q++)
            {
              for (int32_t r = 0; r <= 2; r++)
              {
                m_VoxelUpdateValues->NEIGHBORHOOD[p][q][r] = 0.0;
                m_VoxelUpdateValues->BOUNDARYFLAG[p][q][r] = 0;
              }
            }
          }
#ifndef CIRCULAR_BOUNDARY_CONDITION

          //For a given (i,j,k) store its 26 point neighborhood
          for (int32_t p = -1; p <= 1; p++)
          {
            for (int32_t q = -1; q <= 1; q++)
            {
              for (int32_t r = -1; r <= 1; r++)
              {
                if(i + p >= 0 && i + p < geometry->N_y) if(j_new + q >= 0 && j_new + q < geometry->N_z) if(k_new + r >= 0 && k_new + r < geometry->N_x)
                {
                  m_VoxelUpdateValues->NEIGHBORHOOD[p + 1][q + 1][r + 1] = geometry->Object->d[q + j_new][r + k_new][p + i];
                  m_VoxelUpdateValues->BOUNDARYFLAG[p + 1][q + 1][r + 1] = 1;
                }
                else
                {
                  m_VoxelUpdateValues->BOUNDARYFLAG[p + 1][q + 1][r + 1] = 0;
                }

              }
            }
          }
#else
          for(p = -1; p <=1; p++)
          for(q = -1; q <= 1; q++)
          for(r = -1; r <= 1;r++)
          {
            tempindex_x = mod(r+k_new,geometry->N_x);
            tempindex_y =mod(p+i,geometry->N_y);
            tempindex_z = mod(q+j_new,geometry->N_z);
            NEIGHBORHOOD[p+1][q+1][r+1] = geometry->Object->d[tempindex_z][tempindex_x][tempindex_y];
            BOUNDARYFLAG[p+1][q+1][r+1]=1;
          }

#endif//circular boundary condition check
          m_VoxelUpdateValues->NEIGHBORHOOD[1][1][1] = 0.0;

#ifdef DEBUG
          if(i == 0 && j == 31 && k == 31)
          {
            printf("***************************\n");
            printf("Geom %lf\n", geometry->Object->d[i][31][31]);
            for (int p = 0; p <= 2; p++)
            {
              for (int q = 0; q <= 2; q++)
              {
                for (int32_t r = 0; r <= 2; r++)
                {
                  printf("%lf\n", m_VoxelUpdateValues->NEIGHBORHOOD[p][q][r]);
                }
              }
            }
          }
#endif
          //Compute theta1 and theta2

          V = geometry->Object->d[j_new][k_new][i]; //Store the present value of the voxel
          THETA1 = 0.0;
          THETA2 = 0.0;

          /*    y = ((DATA_TYPE)i+0.5)*Geometry->delta_xy + Geometry->y0;
           t = y;
           tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
           tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax? t + Geometry->delta_xy/2 : Sinogram->TMax;

           slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
           slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

           if(slice_index_min < 0)
           slice_index_min = 0;
           if(slice_index_max >= Sinogram->N_t)
           slice_index_max = Sinogram->N_t-1;*/

          //TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
          for (uint32_t q = 0; q < TempMemBlock->count; q++)
          {

            uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (sinogram->N_r)));
            uint16_t i_r = (TempMemBlock->index[q] % (sinogram->N_r));
            VoxelLineAccessCounter = 0;
            for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
            // for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
            {
              /* center_t = ((DATA_TYPE)i_t + 0.5)*Sinogram->delta_t + Sinogram->T0;
               delta_t = fabs(center_t - t);
               index_delta_t = floor(delta_t/OffsetT);

               if(index_delta_t < DETECTOR_RESPONSE_BINS)
               {
               w3 = delta_t - index_delta_t*OffsetT;
               w4 = (index_delta_t+1)*OffsetT - delta_t;
               //TODO: interpolation
               ProfileThickness =(w4/OffsetT)*H_t[0][i_theta][index_delta_t] + (w3/OffsetT)*H_t[0][i_theta][index_delta_t+1 < DETECTOR_RESPONSE_BINS ? index_delta_t+1:DETECTOR_RESPONSE_BINS-1];
               }
               else
               {
               ProfileThickness=0;
               }*/

              THETA2 += ((NuisanceParams->I_0->d[i_theta] * NuisanceParams->I_0->d[i_theta])
                  * (VoxelLineResponse[i].values[VoxelLineAccessCounter] * VoxelLineResponse[i].values[VoxelLineAccessCounter]) * (TempMemBlock->values[q])
                  * (TempMemBlock->values[q]) * Weight->d[i_theta][i_r][i_t]);
              THETA1 += NuisanceParams->I_0->d[i_theta] * ErrorSino->d[i_theta][i_r][i_t] * (TempMemBlock->values[q])
                  * (VoxelLineResponse[i].values[VoxelLineAccessCounter]) * Weight->d[i_theta][i_r][i_t];
              VoxelLineAccessCounter++;
            }
          }
          THETA1 *= -1;
          minMax(&low, &high);

#ifdef DEBUG
          if(i == 0 && j == 31 && k == 31) printf("(%lf,%lf,%lf) \n", low, high, V - (THETA1 / THETA2));
#endif

          //Solve the 1-D optimization problem
          //printf("V before updating %lf",V);
#ifndef SURROGATE_FUNCTION
          //TODO : What if theta1 = 0 ? Then this will give error
          DerivOfCostFunc docf(m_VoxelUpdateValues->BOUNDARYFLAG, m_VoxelUpdateValues->NEIGHBORHOOD,
                               m_VoxelUpdateValues->FILTER, V, THETA1, THETA2, SIGMA_X_P, MRF_P);
          DATA_TYPE UpdatedVoxelValue;
          UpdatedVoxelValue = (DATA_TYPE)solve < DerivOfCostFunc > (&docf, (double)low, (double)high, (double)accuracy, &errorcode);

#else

          errorcode=0;
#ifdef QGGMRF
          UpdatedVoxelValue = CE_FunctionalSubstitution(low,high);

#else
          SurrogateUpdate=surrogateFunctionBasedMin();
          UpdatedVoxelValue=SurrogateUpdate;
#endif //QGGMRF
#endif//Surrogate function
          //printf("%lf\n",SurrogateUpdate);

          if(errorcode == 0)
          {
            //    printf("(%lf,%lf,%lf)\n",low,high,UpdatedVoxelValue);
            //  printf("Updated %lf\n",UpdatedVoxelValue);
#ifdef POSITIVITY_CONSTRAINT
            if(UpdatedVoxelValue < 0.0)
            { //Enforcing positivity constraints
              UpdatedVoxelValue = 0.0;}
#endif
          }
          else
          {
            if(THETA1 == 0 && low == 0 && high == 0)
            {
              UpdatedVoxelValue = 0;
            }
            else
            {
              printf("Voxel Update Error %d %d\n", j_new, k_new);
            }
          }

          //TODO Print appropriate error messages for other values of error code
          geometry->Object->d[j_new][k_new][i] = UpdatedVoxelValue;

#ifdef ROI
          if(Mask->d[j_new][k_new] == 1)
          {

            AverageUpdate += fabs(geometry->Object->d[j_new][k_new][i] - V);
            AverageMagnitudeOfRecon += fabs(V); //computing the percentage update =(Change in mag/Initial magnitude)
          }
#endif

          //Update the ErrorSinogram

          for (uint32_t q = 0; q < TempMemBlock->count; q++)
          {

            uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (sinogram->N_r)));
            uint16_t i_r = (TempMemBlock->index[q] % (sinogram->N_r));
            VoxelLineAccessCounter = 0;
            for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
            //for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
            {
              /*  center_t = ((DATA_TYPE)i_t + 0.5)*Sinogram->delta_t + Sinogram->T0;
               delta_t = fabs(center_t - t);
               index_delta_t = floor(delta_t/OffsetT);

               if(index_delta_t < DETECTOR_RESPONSE_BINS)
               {
               w3 = delta_t - index_delta_t*OffsetT;
               w4 = (index_delta_t+1)*OffsetT - delta_t;
               //TODO: interpolation
               ProfileThickness =(w4/OffsetT)*H_t[0][i_theta][index_delta_t] + (w3/OffsetT)*H_t[0][i_theta][index_delta_t+1 < DETECTOR_RESPONSE_BINS ? index_delta_t+1:DETECTOR_RESPONSE_BINS-1];
               }
               else
               {
               ProfileThickness=0;
               }*/

              ErrorSino->d[i_theta][i_r][i_t] -= (NuisanceParams->I_0->d[i_theta]
                  * (TempMemBlock->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (geometry->Object->d[j_new][k_new][i] - V)));
              VoxelLineAccessCounter++;
            }
          }
          Idx++;
        }

      }
      else
      {
        continue;
      }

    }
  }
  STOP;
  PRINTTIME;

#ifdef RANDOM_ORDER_UPDATES
  for (int32_t j = 0; j < geometry->N_z; j++)
  { //Row index
    for (int32_t k = 0; k < geometry->N_x; k++)
    { //Column index
      if(VisitCount->d[j][k] == 0)
      {
        printf("Pixel (%d %d) not visited\n", j, k);
      }
    }
  }
#endif

#ifdef COST_CALCULATE
  /*********************Cost Calculation***************************************************/
  cost->d[*cost_counter] = computeCost(ErrorSino, Weight);
  printf("%lf\n", cost->d[*cost_counter]);

  if(cost->d[*cost_counter] - cost->d[*cost_counter - 1] > 0)
  {
    printf("Cost just increased!\n");
    return;
  }

  fwrite(&cost->d[*cost_counter], sizeof(DATA_TYPE), 1, Fp2);
  ++(*cost_counter);
  /*******************************************************************************/
#else
  printf("%d\n", Iter);
#endif //Cost calculation endif
#ifdef ROI
  if(AverageMagnitudeOfRecon > 0)
  {
    printf("%d,%lf\n", Iter + 1, AverageUpdate / AverageMagnitudeOfRecon);
    if((AverageUpdate / AverageMagnitudeOfRecon) < getInputs()->StopThreshold)
    {
      printf("This is the terminating point %d\n", Iter);
      getInputs()->StopThreshold*=THRESHOLD_REDUCTION_FACTOR; //Reducing the thresold for subsequent iterations
      return;
    }
  }
#endif//ROI end


#ifdef WRITE_INTERMEDIATE_RESULTS

  if(Iter == NumOfWrites*WriteCount)
  {
    WriteCount++;
    sprintf(buffer,"%d",Iter);
    sprintf(Filename,"ReconstructedObjectAfterIter");
    strcat(Filename,buffer);
    strcat(Filename,".bin");
    Fp3 = fopen(Filename, "w");
    //  for (i=0; i < Geometry->N_y; i++)
    //    for (j=0; j < Geometry->N_z; j++)
    //      for (k=0; k < Geometry->N_x; k++)
    TempPointer = geometry->Object;
    NumOfBytesWritten=fwrite(&(geometry->Object->d[0][0][0]), sizeof(DATA_TYPE),geometry->N_x*geometry->N_y*geometry->N_z, Fp3);
    printf("%d\n",NumOfBytesWritten);

    fclose(Fp3);
  }
#endif


}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DATA_TYPE VoxelUpdate::computeCost(RealVolumeType::Pointer ErrorSino,RealVolumeType::Pointer Weight)
{
  DATA_TYPE cost=0,temp=0;

  Sinogram* sinogram = getSinogram();
  Geometry* geometry = getGeometry();

  DATA_TYPE MRF_P = m_VoxelUpdateValues->MRF_P;

//Data Mismatch Error
  for (uint16_t i = 0; i < sinogram->N_theta; i++)
  {
    for (uint16_t j = 0; j < sinogram->N_r; j++)
    {
      for (uint16_t k = 0; k < sinogram->N_t; k++)
      {
        cost += (ErrorSino->d[i][j][k] * ErrorSino->d[i][j][k] * Weight->d[i][j][k]);
      }
    }
  }

  cost /= 2;

  printf("Data mismatch term =%lf\n",cost);

//Prior Model Error
#ifndef QGGMRF
  for (uint16_t i = 0; i < geometry->N_z; i++)
    for (uint16_t j = 0; j < geometry->N_x; j++)
      for(uint16_t k = 0; k < geometry->N_y; k++)
      {

        if(k+1 <  geometry->N_y)
          temp += m_VoxelUpdateValues->FILTER[2][1][1]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i][j][k+1]),MRF_P);


        if(j+1 < geometry->N_x)
        {
          if(k-1 >= 0)
            temp += m_VoxelUpdateValues->FILTER[0][1][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k-1]),MRF_P);


          temp += m_VoxelUpdateValues->FILTER[1][1][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k]),MRF_P);


          if(k+1 < geometry->N_y)
            temp += m_VoxelUpdateValues->FILTER[2][1][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k+1]),MRF_P);

        }

        if(i+1 < geometry->N_z)
        {

          if(j-1 >= 0)
            temp += m_VoxelUpdateValues->FILTER[1][2][0]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k]),MRF_P);

          temp += m_VoxelUpdateValues->FILTER[1][2][1]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k]),MRF_P);

          if(j+1 < geometry->N_x)
            temp += m_VoxelUpdateValues->FILTER[1][2][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k]),MRF_P);


          if(j-1 >= 0)
          {
            if(k-1 >= 0)
              temp += m_VoxelUpdateValues->FILTER[0][2][0]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k-1]),MRF_P);

            if(k+1 < geometry->N_y)
              temp += m_VoxelUpdateValues->FILTER[2][2][0]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k+1]),MRF_P);

          }

          if(k-1 >= 0)
            temp += m_VoxelUpdateValues->FILTER[0][2][1]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k-1]),MRF_P);

          if(j+1 < geometry->N_x)
          {
            if(k-1 >= 0)
              temp += m_VoxelUpdateValues->FILTER[0][2][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k-1]),MRF_P);

            if(k+1 < geometry->N_y)
              temp+= m_VoxelUpdateValues->FILTER[2][2][2]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k+1]),MRF_P);
          }

          if(k+1 < geometry->N_y)
            temp+= m_VoxelUpdateValues->FILTER[2][2][1]*pow(fabs(geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k+1]),MRF_P);
        }
      }
    cost+=(temp/(m_VoxelUpdateValues->MRF_P*m_VoxelUpdateValues->SIGMA_X_P));
#else
  for (i = 0; i < geometry->N_z; i++)
    for (j = 0; j < geometry->N_x; j++)
      for(k = 0; k < geometry->N_y; k++)
      {

        if(k+1 <  geometry->N_y)
        {
          delta=geometry->Object->d[i][j][k]-geometry->Object->d[i][j][k+1];
          temp += FILTER[2][1][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
        }



        if(j+1 < geometry->N_x)
        {
          if(k-1 >= 0)
          {
            delta=geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k-1];
            temp += FILTER[0][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }

          delta=geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k];
          temp += FILTER[1][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));


          if(k+1 < geometry->N_y)
          {
            delta=geometry->Object->d[i][j][k]-geometry->Object->d[i][j+1][k+1];
            temp += FILTER[2][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }

        }

        if(i+1 < geometry->N_z)
        {

          if(j-1 >= 0)
          {
            delta = geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k];
            temp += FILTER[1][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }

          delta=geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k];
          temp += FILTER[1][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));

          if(j+1 < geometry->N_x)
          {
            delta=geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k];
            temp += FILTER[1][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }


          if(j-1 >= 0)
          {
            if(k-1 >= 0)
            {
              delta=geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k-1];
              temp += FILTER[0][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
            }

            if(k+1 < geometry->N_y)
            {
              delta=geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j-1][k+1];
              temp += FILTER[2][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
            }

          }

          if(k-1 >= 0)
          {
            delta = geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k-1];
            temp += FILTER[0][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }

          if(j+1 < geometry->N_x)
          {
            if(k-1 >= 0)
            {
              delta = geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k-1];
              temp += FILTER[0][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
            }

            if(k+1 < geometry->N_y)
            {
              delta = geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j+1][k+1];
              temp+= FILTER[2][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
            }
          }

          if(k+1 < geometry->N_y)
          {
            delta = geometry->Object->d[i][j][k]-geometry->Object->d[i+1][j][k+1];
            temp+= FILTER[2][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
          }
        }
      }
      cost+=(temp);
#endif //QGGMRF

  //printf("Cost calculation End..\n");


//Noise Error
#ifdef NOISE_MODEL
  temp=0;
  for(uint16_t i=0;i< getSinogram()->N_theta;i++) {
    if(Weight->d[i][0][0] != 0) {
    temp += log(2*M_PI*(1.0/Weight->d[i][0][0]));//2*pi*sigma_k^{2}
    }}

  temp*=((getSinogram()->N_r * getSinogram()->N_t)/2);

  cost+=temp;
#endif//noise model
  return cost;
}
