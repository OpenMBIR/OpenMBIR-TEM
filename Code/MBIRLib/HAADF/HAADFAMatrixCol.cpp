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


#include "HAADFAMatrixCol.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFAMatrixCol::HAADFAMatrixCol(size_t* dims, int32_t c)
{
  valuesPtr = RealArrayType::New(dims, "VoxelLineResponse_Values");
  values = valuesPtr->getPointer(0);
  indexPtr = UInt32ArrayType::New(dims, "VoxelLineResponse_index");
  index = indexPtr->getPointer(0);
  count = c;
  d0 = 0xABABABABABABABABull;
  d1 = 0xCACACACACACACACAull;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFAMatrixCol::~HAADFAMatrixCol()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFAMatrixCol::setCount(uint32_t c)
{
  if(c > valuesPtr->getDims()[0])
  {
    std::cout << "BAD!!! c: " << c << "  count: " << count << std::endl;
    assert(false);
  }
  count = c;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFAMatrixCol::Pointer HAADFAMatrixCol::calculateHAADFAMatrixColumnPartial(SinogramPtr sinogram,
                                                                             GeometryPtr geometry,
                                                                             TomoInputsPtr tomoInputs,
                                                                             AdvancedParametersPtr advParams,
                                                                             uint16_t row,
                                                                             uint16_t col,
                                                                             uint16_t slice,
                                                                             RealVolumeType::Pointer detectorResponse,
                                                                             HAADFDetectorParameters::Pointer haadfParameters)
{


  int32_t j, k, sliceidx;
  Real_t x, z, y;
  Real_t r; //this is used to find where does the ray passing through the voxel at certain angle hit the detector
  Real_t t; //this is similar to r but along the y direction
  Real_t tmin, tmax;
  Real_t rmax, rmin; //stores the start and end points of the pixel profile on the detector
  Real_t R_Center, TempConst, checksum = 0, delta_r;
//  DATA_TYPE Integral = 0;
  Real_t T_Center, delta_t;
  Real_t MaximumSpacePerColumn; //we will use this to allocate space
  Real_t AvgNumXElements, AvgNumYElements; //This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
//  DATA_TYPE ProfileThickness,stepsize;

  //interpolation variables
  Real_t w1, w2, w3, w4, f1, InterpolatedValue, ContributionAlongT;
//  DATA_TYPE f2;
  int32_t index_min, index_max, slice_index_min, slice_index_max, index_delta_r, index_delta_t; //stores the detector index in which the profile lies
  int32_t BaseIndex, FinalIndex;
//  int32_t ProfileIndex=0;
//  int32_t NumOfDisplacements=32;
  uint32_t count = 0;

  sliceidx = 0;



  x = geometry->x0 + ((Real_t)col + 0.5) * tomoInputs->delta_xz; //0.5 is for center of voxel. x_0 is the left corner
  z = geometry->z0 + ((Real_t)row + 0.5) * tomoInputs->delta_xz; //0.5 is for center of voxel. x_0 is the left corner
  y = geometry->y0 + ((Real_t)slice + 0.5) * tomoInputs->delta_xy;

  TempConst = (advParams->PROFILE_RESOLUTION) / (2 * tomoInputs->delta_xz);

  //alternately over estimate the maximum size require for a single AMatrix column
  AvgNumXElements = ceil(3 * tomoInputs->delta_xz / sinogram->delta_r);
  AvgNumYElements = ceil(3 * tomoInputs->delta_xy / sinogram->delta_t);
  MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements) * sinogram->N_theta;

  size_t dims[1] = { MaximumSpacePerColumn };
  HAADFAMatrixCol::Pointer Temp = HAADFAMatrixCol::New(dims, 0);
//  HAADFAMatrixCol* Temp = (HAADFAMatrixCol*)get_spc(1, sizeof(HAADFAMatrixCol)); //This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end
//
//  Temp->values = (Real_t*)get_spc((uint32_t)MaximumSpacePerColumn, sizeof(Real_t));
//  Temp->index = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn, sizeof(uint32_t));

  RealArrayType::Pointer cosine = haadfParameters->getcosine();
  RealArrayType::Pointer sine = haadfParameters->getsine();
  const Real_t OffsetR = haadfParameters->getOffsetR();
  const Real_t OffsetT = haadfParameters->getOffsetT();

  if(advParams->AREA_WEIGHTED)
  {
    for (uint32_t i = 0; i < sinogram->N_theta; i++)
    {

      r = x * cosine->d[i] - z * sine->d[i];
      t = y;

      rmin = r - tomoInputs->delta_xz;
      rmax = r + tomoInputs->delta_xz;

      tmin = (t - tomoInputs->delta_xy / 2) > sinogram->T0 ? t - tomoInputs->delta_xy / 2 : sinogram->T0;
      tmax = (t + tomoInputs->delta_xy / 2) <= sinogram->TMax ? t + tomoInputs->delta_xy / 2 : sinogram->TMax;

      if(rmax < sinogram->R0 || rmin > sinogram->RMax) continue;

      index_min = static_cast<int32_t>(floor(((rmin - sinogram->R0) / sinogram->delta_r)));
      index_max = static_cast<int32_t>(floor((rmax - sinogram->R0) / sinogram->delta_r));

      if(index_max >= sinogram->N_r) index_max = sinogram->N_r - 1;

      if(index_min < 0) index_min = 0;

      slice_index_min = static_cast<int32_t>(floor((tmin - sinogram->T0) / sinogram->delta_t));
      slice_index_max = static_cast<int32_t>(floor((tmax - sinogram->T0) / sinogram->delta_t));

      if(slice_index_min < 0) slice_index_min = 0;
      if(slice_index_max >= sinogram->N_t) slice_index_max = sinogram->N_t - 1;

      BaseIndex = i * sinogram->N_r; //*Sinogram->N_t;

      for (j = index_min; j <= index_max; j++) //Check
      {

        //Accounting for Beam width
        R_Center = (sinogram->R0 + (((Real_t)j) + 0.5) * (sinogram->delta_r)); //the 0.5 is to get to the center of the detector

        //Find the difference between the center of detector and center of projection and compute the Index to look up into
        delta_r = fabs(r - R_Center);
        index_delta_r = static_cast<int32_t>(floor((delta_r / OffsetR)));

        if(index_delta_r >= 0 && index_delta_r < advParams->DETECTOR_RESPONSE_BINS)
        {
          T_Center = (sinogram->T0 + (((Real_t)sliceidx) + 0.5) * (sinogram->delta_t));
          delta_t = fabs(t - T_Center);
          index_delta_t = 0; //floor(delta_t/OffsetT);
          if(index_delta_t >= 0 && index_delta_t < advParams->DETECTOR_RESPONSE_BINS)
          {
            //Using index_delta_t,index_delta_t+1,index_delta_r and index_delta_r+1 do bilinear interpolation
            w1 = delta_r - index_delta_r * OffsetR;
            w2 = (index_delta_r + 1) * OffsetR - delta_r;

            w3 = delta_t - index_delta_t * OffsetT;
            w4 = (index_delta_r + 1) * OffsetT - delta_t;

            uint16_t iidx = index_delta_r + 1 < advParams->DETECTOR_RESPONSE_BINS ? index_delta_r + 1 : advParams->DETECTOR_RESPONSE_BINS - 1;
            f1 = (w2 / OffsetR) * detectorResponse->getValue(index_delta_t, i, index_delta_r)
                + (w1 / OffsetR) * detectorResponse->getValue(index_delta_t, i, iidx);

            if(sliceidx == slice_index_min) ContributionAlongT = (sliceidx + 1) * sinogram->delta_t - tmin;
            else if(sliceidx == slice_index_max) ContributionAlongT = tmax - (sliceidx) * sinogram->delta_t;
            else
            {
              ContributionAlongT = sinogram->delta_t;
            }
            InterpolatedValue = f1; //*ContributionAlongT;//(w3/OffsetT)*f2 + (w4/OffsetT)*f2;
            if(InterpolatedValue > 0)
            {
              FinalIndex = BaseIndex + (int32_t)j; //+ (int32_t)sliceidx * Sinogram->N_r;
              Temp->values[count] = InterpolatedValue; //detectorResponse[index_delta_t][i][index_delta_r];
              Temp->index[count] = FinalIndex; //can instead store a triple (row,col,slice) for the sinogram
              count++;
            }
          }
        }
      }
    }
  }

  dims[0] = count;
  HAADFAMatrixCol::Pointer Ai = HAADFAMatrixCol::New(dims, 0);

  k = 0;
  for (uint32_t i = 0; i < count; i++)
  {
    if(Temp->values[i] > 0.0)
    {
      Ai->values[k] = Temp->values[i];
      checksum += Ai->values[k];
      Ai->index[k] = Temp->index[i];
      k++;
    }
  }
  Ai->setCount(k);

  return Ai;
}
