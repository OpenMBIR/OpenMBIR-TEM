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
#include "CalculateAMatrixColumn.h"

#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CalculateAMatrixColumn::CalculateAMatrixColumn()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CalculateAMatrixColumn::~CalculateAMatrixColumn()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CalculateAMatrixColumn::execute()
{

  int32_t j,k,sliceidx;
  Real_t x,z,y;
  Real_t r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
  Real_t t; //this is similar to r but along the y direction
  Real_t tmin,tmax;
  Real_t rmax,rmin;//stores the start and end points of the pixel profile on the detector
  Real_t RTemp,TempConst,checksum = 0,Integral = 0;
  Real_t LeftEndOfBeam;
  Real_t MaximumSpacePerColumn;//we will use this to allocate space
  Real_t AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
  Real_t ProfileThickness;
  int32_t index_min,index_max,slice_index_min,slice_index_max;//stores the detector index in which the profile lies
  int32_t BaseIndex,FinalIndex,ProfileIndex=0;
  uint32_t count = 0;

  SinogramPtr sinogram = getSinogram();
  GeometryPtr geometry = getGeometry();
  TomoInputsPtr inputs = getTomoInputs();

#ifdef BEAM_CALCULATION
  BEAM_WIDTH = (0.5)*sinogram->delta_r;
#else
  BEAM_WIDTH = sinogram->delta_r;
#endif

#ifdef DISTANCE_DRIVEN
  Real_t d1,d2; //These are the values of the detector boundaries
#endif

  Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
  AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

  // printf("Space allocated for column %d %d\n",row,col);

  //Temp->index = (uint32_t*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(uint32_t));
  //Temp->values = (DATA_TYPE*)multialloc(sizeof(DATA_TYPE),1,Sinogram->N_r*Sinogram->N_theta);//makes the values =0

  x = geometry->x0 + ((Real_t)col+0.5)*inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  z = geometry->z0 + ((Real_t)row+0.5)*inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  y = geometry->y0 + ((Real_t)slice + 0.5)*inputs->delta_xy;

  TempConst=(PROFILE_RESOLUTION)/(2*inputs->delta_xz);


  //  Temp->values = (DATA_TYPE*)calloc(Sinogram->N_t*Sinogram->N_r*Sinogram->N_theta,sizeof(DATA_TYPE));//(DATA_TYPE*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(DATA_TYPE));//makes the values =0

  //alternately over estimate the maximum size require for a single AMatrix column
  AvgNumXElements = ceil(3*inputs->delta_xz/sinogram->delta_r);
  AvgNumYElements = ceil(3*inputs->delta_xy/sinogram->delta_t);
  MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*sinogram->N_theta;

  Temp->values = (Real_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(Real_t));
  Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));

  //printf("%lf",Temp->values[10]);

#ifdef AREA_WEIGHTED
  for(uint32_t i=0;i<sinogram->N_theta;i++)
  {

    r = x*cosine[i] - z*sine[i];
    t = y;

    rmin = r - inputs->delta_xz;
    rmax = r + inputs->delta_xz;

    tmin = (t - inputs->delta_xy/2) > sinogram->T0 ? t-inputs->delta_xy/2 : sinogram->T0;
    tmax = (t + inputs->delta_xy/2) <= sinogram->TMax ? t + inputs->delta_xy/2 : sinogram->TMax;

    if(rmax < sinogram->R0 || rmin > sinogram->RMax)
      continue;



    index_min = floor(((rmin - sinogram->R0)/sinogram->delta_r));
    index_max = floor((rmax - sinogram->R0)/sinogram->delta_r);

    if(index_max >= sinogram->N_r)
      index_max = sinogram->N_r - 1;

    if(index_min < 0)
      index_min = 0;

    slice_index_min = floor((tmin - sinogram->T0)/sinogram->delta_t);
    slice_index_max = floor((tmax - sinogram->T0)/sinogram->delta_t);

    if(slice_index_min < 0)
      slice_index_min = 0;
    if(slice_index_max >= sinogram->N_t)
      slice_index_max = sinogram->N_t -1;

    BaseIndex = i*sinogram->N_r*sinogram->N_t;

    // if(row == 63 && col == 63)
    //    printf("%d %d\n",index_min,index_max);


    //Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
    //of same dimesions


    for(j = index_min;j <= index_max; j++)//Check
    {

      Integral = 0.0;

      //Accounting for Beam width
      RTemp = (sinogram->R0 + (((Real_t)j) + 0.5) *(sinogram->delta_r));//the 0.5 is to get to the center of the detector

      LeftEndOfBeam = RTemp - (BEAM_WIDTH/2);//Consider pre storing the divide by 2

      //   if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
      //printf("%d %lf %lf\n",j, LeftEndOfBeam,rmin);


      for(k=0; k < BEAM_RESOLUTION; k++)
      {

        RTemp = LeftEndOfBeam + ((((Real_t)k)*(BEAM_WIDTH))/BEAM_RESOLUTION);

        if (RTemp-rmin >= 0.0)
        {
          ProfileIndex = (int32_t)floor((RTemp-rmin)*TempConst);//Finding the nearest neighbor profile to the beam
          //if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
          //printf("%d\n",ProfileIndex);
          if(ProfileIndex > PROFILE_RESOLUTION)
            ProfileIndex = PROFILE_RESOLUTION;
        }
        if(ProfileIndex < 0)
          ProfileIndex = 0;


        if(ProfileIndex >= 0 && ProfileIndex < PROFILE_RESOLUTION)
        {
#ifdef BEAM_CALCULATION
          Integral+=(BeamProfile[k]*VoxelProfile[i][ProfileIndex]);
#else
          Integral+=(VoxelProfile[i][ProfileIndex]/PROFILE_RESOLUTION);
#endif
          //  if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
          //   printf("Index %d %lf Voxel %lf I=%d\n",ProfileIndex,BeamProfile[k],VoxelProfile[2][274],i);
        }

      }
      if(Integral > 0.0)
      {
        //  printf("Entering, Final Index %d %d\n",FinalIndex,Temp->values[0]);
        //    printf("Done %d %d\n",slice_index_min,slice_index_max);
        for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
        {
          if(sliceidx == slice_index_min)
            ProfileThickness = (((sliceidx+1)*sinogram->delta_t + sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
          else
          {
            if (sliceidx == slice_index_max)
            {
              ProfileThickness = (tmax - ((sliceidx)*sinogram->delta_t + sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
            }
            else
            {
              ProfileThickness = sinogram->delta_t;//Sinogram->delta_t;
            }

          }
          if(ProfileThickness > 0)
          {

            FinalIndex = BaseIndex + (int32_t)j + (int32_t)sliceidx * sinogram->N_r;
            Temp->values[count] = Integral*ProfileThickness;
            Temp->index[count] = FinalIndex;//can instead store a triple (row,col,slice) for the sinogram
            //printf("Done\n");
#ifdef CORRECTION
            Temp->values[count]/=NORMALIZATION_FACTOR[FinalIndex];//There is a normalizing constant for each measurement . So we have a total of Sinogram.N_x * Sinogram->N_theta values
#endif
            //printf("Access\n");
            count++;
          }

        }
      }

      //  End of Beam width accounting


      //        ProfileIndex=(uint32_t)(((RTemp-rmin)*TempConst));
      // //    //   printf("%d\n",ProfileIndex);
      //        if(ProfileIndex>=0 && ProfileIndex < PROFILE_RESOLUTION)
      //        {
      //    if(VoxelProfile[i][ProfileIndex] > 0.0)
      //    {
      //          Temp->values[FinalIndex]=VoxelProfile[i][ProfileIndex];
      //         // Temp->index[count++]=FinalIndex;
      //   count++;
      //    }
      //        }
    }



  }

#endif

#ifdef DISTANCE_DRIVEN

  for(uint16_t i=0;i<sinogram->N_theta;i++)
  {

    r = x*cosine[i] - z*sine[i];
    t = y;

    tmin = (t - inputs->delta_xy/2) > -geometry->LengthY/2 ? t-inputs->delta_xy/2 : -geometry->LengthY/2;
    tmax = (t + inputs->delta_xy/2) <= geometry->LengthY/2 ? t + inputs->delta_xy/2 : geometry->LengthY/2;


    if(sinogram->angles[i]*(180/M_PI) >= -45 && sinogram->angles[i]*(180/M_PI) <= 45)
    {
      rmin = r - (inputs->delta_xz/2)*(cosine[i]);
      rmax = r + (inputs->delta_xz/2)*(cosine[i]);
    }

    else
    {
      rmin = r - (inputs->delta_xz/2)*fabs(sine[i]);
      rmax = r + (inputs->delta_xz/2)*fabs(sine[i]);
    }


    if(rmax < sinogram->R0 || rmin > sinogram->R0 + sinogram->N_r*sinogram->delta_r)
      continue;

    index_min = floor(((rmin-sinogram->R0)/sinogram->delta_r));
    index_max = floor((rmax-sinogram->R0)/sinogram->delta_r);

    slice_index_min = floor((tmin - sinogram->T0)/sinogram->delta_t);
    slice_index_max = floor((tmax - sinogram->T0)/sinogram->delta_t);

    if(slice_index_min < 0)
      slice_index_min = 0;
    if(slice_index_max >= sinogram->N_t)
      slice_index_max = sinogram->N_t -1;


    if(index_max >= sinogram->N_r)
      index_max = sinogram->N_r-1;

    if(index_min < 0)
      index_min = 0;

    BaseIndex = i*sinogram->N_r*sinogram->N_t;

    // if(row == 63 && col == 63)
    //    printf("%d %d\n",index_min,index_max);



    //Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
    //of same dimesions

    for(j = index_min;j <= index_max; j++)//Check
    {
      d1  = ((Real_t)j)*sinogram->delta_r + sinogram->R0;
      d2 =  d1 + sinogram->delta_r;

      if(rmax < d1)
      {
        Integral = 0;
      }
      else
      {
        if(rmin > d1 && rmin < d2 && rmax > d2)
        {
          Integral = (d2 - rmin)/sinogram->delta_r;
        }
        else
        {
          if(rmin >= d1 && rmin <= d2 && rmax >= d1 && rmax <= d2)
            Integral= (rmax - rmin)/sinogram->delta_r;
          else
            if(rmax > d1 && rmax < d2 && rmin < d1)
              Integral = (rmax - d1)/sinogram->delta_r;
            else
              if( rmin < d1 && rmax > d2)
                Integral = 1;
              else
                Integral = 0;
        }


      }
      if(Integral > 0)
      {
        //   printf("Final Index %d %lf\n",FinalIndex,cosine[i]);
        for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
        {
          if(sliceidx == slice_index_min)
            ProfileThickness = (((sliceidx+1)*sinogram->delta_t + sinogram->T0) - tmin)*sinogram->delta_t;
          else {
            if (sliceidx == slice_index_max)
            {
              ProfileThickness = (tmax - ((sliceidx)*sinogram->delta_t + sinogram->T0))*sinogram->delta_t;
            }
            else {
              ProfileThickness = sinogram->delta_t;
            }

          }
          if (ProfileThickness > 0)
          {
            FinalIndex = BaseIndex + (uint32_t)j + (int32_t)sliceidx * sinogram->N_r;

            Temp->values[count] = Integral*ProfileThickness;
            Temp->index[count] = FinalIndex;
#ifdef CORRECTION
            Temp->values[count]/=NORMALIZATION_FACTOR[FinalIndex];//There is a normalizing constant for each measurement . So we have a total of Sinogram.N_x * Sinogram->N_theta values
#endif
            count++;
          }
        }
      }
      else
      {
        //  printf("%lf \n",Sinogram->angles[i]*180/PI);
      }
    }



  }


#endif
  // printf("Final Space allocation for column %d %d\n",row,col);

  Ai->values=(Real_t*)get_spc(count,sizeof(Real_t));
  Ai->index=(uint32_t*)get_spc(count,sizeof(uint32_t));
  k=0;
  for(uint32_t i = 0; i < count; i++)
  {
    if(Temp->values[i] > 0.0)
    {
      Ai->values[k]=Temp->values[i];
      checksum+=Ai->values[k];
      Ai->index[k++]=Temp->index[i];
    }

  }

  Ai->count=k;

  //printf("%d %d \n",Ai->count,count);

  //printf("(%d,%d) %lf \n",row,col,checksum);

  free(Temp->values);
  free(Temp->index);
  free(Temp);

}
