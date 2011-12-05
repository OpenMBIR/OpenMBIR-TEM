/*
 * CalculateAMatrix.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "CalculateAMatrixColumn.h"

#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CalculateAMatrixColumn::CalculateAMatrixColumn() :
m_Inputs(NULL),
m_Sinogram(NULL)
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
  DATA_TYPE x,z,y;
  DATA_TYPE r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
  DATA_TYPE t; //this is similar to r but along the y direction
  DATA_TYPE tmin,tmax;
  DATA_TYPE rmax,rmin;//stores the start and end points of the pixel profile on the detector
  DATA_TYPE RTemp,TempConst,checksum = 0,Integral = 0;
  DATA_TYPE LeftEndOfBeam;
  DATA_TYPE MaximumSpacePerColumn;//we will use this to allocate space
  DATA_TYPE AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
  DATA_TYPE ProfileThickness;
  int32_t index_min,index_max,slice_index_min,slice_index_max;//stores the detector index in which the profile lies
  int32_t BaseIndex,FinalIndex,ProfileIndex=0;
  uint32_t count = 0;


#ifdef BEAM_CALCULATION
  BEAM_WIDTH = (0.5)*m_Sinogram->delta_r;
#else
  BEAM_WIDTH = m_Sinogram->delta_r;
#endif

#ifdef DISTANCE_DRIVEN
  DATA_TYPE d1,d2; //These are the values of the detector boundaries
#endif

  Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
  AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

  // printf("Space allocated for column %d %d\n",row,col);

  //Temp->index = (uint32_t*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(uint32_t));
  //Temp->values = (DATA_TYPE*)multialloc(sizeof(DATA_TYPE),1,Sinogram->N_r*Sinogram->N_theta);//makes the values =0

  x = m_Geometry->x0 + ((DATA_TYPE)col+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  z = m_Geometry->z0 + ((DATA_TYPE)row+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  y = m_Geometry->y0 + ((DATA_TYPE)slice + 0.5)*m_Inputs->delta_xy;

  TempConst=(PROFILE_RESOLUTION)/(2*m_Inputs->delta_xz);


  //  Temp->values = (DATA_TYPE*)calloc(Sinogram->N_t*Sinogram->N_r*Sinogram->N_theta,sizeof(DATA_TYPE));//(DATA_TYPE*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(DATA_TYPE));//makes the values =0

  //alternately over estimate the maximum size require for a single AMatrix column
  AvgNumXElements = ceil(3*m_Inputs->delta_xz/m_Sinogram->delta_r);
  AvgNumYElements = ceil(3*m_Inputs->delta_xy/m_Sinogram->delta_t);
  MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*m_Sinogram->N_theta;

  Temp->values = (DATA_TYPE*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(DATA_TYPE));
  Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));

  //printf("%lf",Temp->values[10]);

#ifdef AREA_WEIGHTED
  for(uint32_t i=0;i<m_Sinogram->N_theta;i++)
  {

    r = x*cosine[i] - z*sine[i];
    t = y;

    rmin = r - m_Inputs->delta_xz;
    rmax = r + m_Inputs->delta_xz;

    tmin = (t - m_Inputs->delta_xy/2) > m_Sinogram->T0 ? t-m_Inputs->delta_xy/2 : m_Sinogram->T0;
    tmax = (t + m_Inputs->delta_xy/2) <= m_Sinogram->TMax ? t + m_Inputs->delta_xy/2 : m_Sinogram->TMax;

    if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax)
      continue;



    index_min = floor(((rmin - m_Sinogram->R0)/m_Sinogram->delta_r));
    index_max = floor((rmax - m_Sinogram->R0)/m_Sinogram->delta_r);

    if(index_max >= m_Sinogram->N_r)
      index_max = m_Sinogram->N_r - 1;

    if(index_min < 0)
      index_min = 0;

    slice_index_min = floor((tmin - m_Sinogram->T0)/m_Sinogram->delta_t);
    slice_index_max = floor((tmax - m_Sinogram->T0)/m_Sinogram->delta_t);

    if(slice_index_min < 0)
      slice_index_min = 0;
    if(slice_index_max >= m_Sinogram->N_t)
      slice_index_max = m_Sinogram->N_t -1;

    BaseIndex = i*m_Sinogram->N_r*m_Sinogram->N_t;

    // if(row == 63 && col == 63)
    //    printf("%d %d\n",index_min,index_max);


    //Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
    //of same dimesions


    for(j = index_min;j <= index_max; j++)//Check
    {

      Integral = 0.0;

      //Accounting for Beam width
      RTemp = (m_Sinogram->R0 + (((DATA_TYPE)j) + 0.5) *(m_Sinogram->delta_r));//the 0.5 is to get to the center of the detector

      LeftEndOfBeam = RTemp - (BEAM_WIDTH/2);//Consider pre storing the divide by 2

      //   if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
      //printf("%d %lf %lf\n",j, LeftEndOfBeam,rmin);


      for(k=0; k < BEAM_RESOLUTION; k++)
      {

        RTemp = LeftEndOfBeam + ((((DATA_TYPE)k)*(BEAM_WIDTH))/BEAM_RESOLUTION);

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
            ProfileThickness = (((sliceidx+1)*m_Sinogram->delta_t + m_Sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
          else
          {
            if (sliceidx == slice_index_max)
            {
              ProfileThickness = (tmax - ((sliceidx)*m_Sinogram->delta_t + m_Sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
            }
            else
            {
              ProfileThickness = m_Sinogram->delta_t;//Sinogram->delta_t;
            }

          }
          if(ProfileThickness > 0)
          {

            FinalIndex = BaseIndex + (int32_t)j + (int32_t)sliceidx * m_Sinogram->N_r;
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

  for(uint16_t i=0;i<m_Sinogram->N_theta;i++)
  {

    r = x*cosine[i] - z*sine[i];
    t = y;

    tmin = (t - m_Inputs->delta_xy/2) > -m_Geometry->LengthY/2 ? t-m_Inputs->delta_xy/2 : -m_Geometry->LengthY/2;
    tmax = (t + m_Inputs->delta_xy/2) <= m_Geometry->LengthY/2 ? t + m_Inputs->delta_xy/2 : m_Geometry->LengthY/2;


    if(m_Sinogram->angles[i]*(180/M_PI) >= -45 && m_Sinogram->angles[i]*(180/M_PI) <= 45)
    {
      rmin = r - (m_Inputs->delta_xz/2)*(cosine[i]);
      rmax = r + (m_Inputs->delta_xz/2)*(cosine[i]);
    }

    else
    {
      rmin = r - (m_Inputs->delta_xz/2)*fabs(sine[i]);
      rmax = r + (m_Inputs->delta_xz/2)*fabs(sine[i]);
    }


    if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->R0 + m_Sinogram->N_r*m_Sinogram->delta_r)
      continue;

    index_min = floor(((rmin-m_Sinogram->R0)/m_Sinogram->delta_r));
    index_max = floor((rmax-m_Sinogram->R0)/m_Sinogram->delta_r);

    slice_index_min = floor((tmin - m_Sinogram->T0)/m_Sinogram->delta_t);
    slice_index_max = floor((tmax - m_Sinogram->T0)/m_Sinogram->delta_t);

    if(slice_index_min < 0)
      slice_index_min = 0;
    if(slice_index_max >= m_Sinogram->N_t)
      slice_index_max = m_Sinogram->N_t -1;


    if(index_max >= m_Sinogram->N_r)
      index_max = m_Sinogram->N_r-1;

    if(index_min < 0)
      index_min = 0;

    BaseIndex = i*m_Sinogram->N_r*m_Sinogram->N_t;

    // if(row == 63 && col == 63)
    //    printf("%d %d\n",index_min,index_max);



    //Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
    //of same dimesions

    for(j = index_min;j <= index_max; j++)//Check
    {
      d1  = ((DATA_TYPE)j)*m_Sinogram->delta_r + m_Sinogram->R0;
      d2 =  d1 + m_Sinogram->delta_r;

      if(rmax < d1)
      {
        Integral = 0;
      }
      else
      {
        if(rmin > d1 && rmin < d2 && rmax > d2)
        {
          Integral = (d2 - rmin)/m_Sinogram->delta_r;
        }
        else
        {
          if(rmin >= d1 && rmin <= d2 && rmax >= d1 && rmax <= d2)
            Integral= (rmax - rmin)/m_Sinogram->delta_r;
          else
            if(rmax > d1 && rmax < d2 && rmin < d1)
              Integral = (rmax - d1)/m_Sinogram->delta_r;
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
            ProfileThickness = (((sliceidx+1)*m_Sinogram->delta_t + m_Sinogram->T0) - tmin)*m_Sinogram->delta_t;
          else {
            if (sliceidx == slice_index_max)
            {
              ProfileThickness = (tmax - ((sliceidx)*m_Sinogram->delta_t + m_Sinogram->T0))*m_Sinogram->delta_t;
            }
            else {
              ProfileThickness = m_Sinogram->delta_t;
            }

          }
          if (ProfileThickness > 0)
          {
            FinalIndex = BaseIndex + (uint32_t)j + (int32_t)sliceidx * m_Sinogram->N_r;

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

  Ai->values=(DATA_TYPE*)get_spc(count,sizeof(DATA_TYPE));
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
