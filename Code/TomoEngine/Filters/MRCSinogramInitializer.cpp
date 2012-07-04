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

#include "MRCSinogramInitializer.h"

#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCSinogramInitializer::MRCSinogramInitializer()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCSinogramInitializer::~MRCSinogramInitializer()
{
}




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCSinogramInitializer::execute()
{
  std::stringstream ss;
  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr inputs = getTomoInputs();
 // int16_t i,j,k;
 // uint16_t TotalNumMaskedViews;

  Real_t sum=0;

  MRCReader::Pointer reader = MRCReader::New(true);
  MRCHeader header;
  int err = reader->readHeader(inputs->sinoFile, &header);
  //reader->printHeader(&header, std::cout);
	if (err < 0)
  {
	  FREE_FEI_HEADERS( header.feiHeaders )
  }

  if (header.mode != 1)
  {
    FREE_FEI_HEADERS( header.feiHeaders )
    ss << "16 bit integers are only supported. Error at line  " << __LINE__ << " in file " << __FILE__ << std::endl;
    setErrorCondition(-1);
    notify(ss.str(), 0, Observable::UpdateErrorMessage);
    return;
  }

  int voxelMin[3] = {0, 0, 0};
  int voxelMax[3] = {header.nx-1, header.ny-1, header.nz-1};
  inputs->fileXSize = header.nx;
  inputs->fileYSize = header.ny;
  inputs->fileZSize = header.nz;

  Real_t CenterOfRot = header.nx/2;
  ss.str("");
  ss << "Center of rotation in this data set is " << CenterOfRot << std::endl;

  if (inputs->useSubvolume == true)
  {
    voxelMin[0] = inputs->xStart;
    voxelMin[1] = inputs->yStart;
    voxelMin[2] = inputs->zStart;

    voxelMax[0] = inputs->xEnd;
    voxelMax[1] = inputs->yEnd;
    voxelMax[2] = inputs->zEnd;

    //Code to ensure region selected has the "right" dimensions
    /************************************************************/
    Real_t LeftLength = CenterOfRot - voxelMin[0];
    Real_t RightLength = voxelMax[0] - CenterOfRot + 1; //1 is present to account for indexing of subvolume which starts from 0

    if(LeftLength < 0)
    {
      voxelMin[0] = CenterOfRot + LeftLength;
      inputs->xStart = voxelMin[0];
      LeftLength *= -1;
    }
    if(RightLength < 0)
    {
      voxelMax[0] = CenterOfRot - RightLength;
      inputs->xEnd = voxelMax[0];
      RightLength *= -1;
    }

    if(LeftLength != RightLength)
    {
      ss << "The subvolume is not symmetric about the center. Adjusting.." << std::endl;
      if(LeftLength < RightLength)
      {
        Real_t tempx = CenterOfRot + LeftLength - 1;
        //Make sure the adjustment does not overrun the size of the data
        if(tempx >= header.nx) tempx = header.nx - 1;

        voxelMax[0] = tempx;
        inputs->xEnd = voxelMax[0];
        ss << "New xEnd : " << voxelMax[0] << std::endl;
      }
      else
      {
        Real_t tempx = CenterOfRot - RightLength;
        //Make sure the adjustment does not overrun the size of the data
        if(tempx < 0) tempx = 0;

        voxelMin[0] = tempx;
        inputs->xStart = voxelMin[0];
        ss << "New xStart : " << voxelMin[0] << std::endl;
      }
    }

    //This part of selecting y can be ignored if the user has selected
    //single slice mode

    //Adjusting the volume along the y-directions so we dont have
    //  issues with pixelation
    int16_t disty = inputs->yEnd - inputs->yStart + 1;
    ss << "Interpolate Factor=" << inputs->interpolateFactor << std::endl;
    //3*iterpFactor is to account for the prior which operates on
    //26 point 3-D neighborhood which needs 3 x-z slices at the least
    int16_t rem_temp = disty % ((int16_t)inputs->interpolateFactor * 3);
    if(rem_temp != 0)
    {
      ss << "The number of y-pixels is not a proper multiple for multi-res" << std::endl;
      int16_t remainder = static_cast<int16_t>((inputs->interpolateFactor * 3) - (rem_temp));

      //Make sure the adjustment does not overrun the size of the data
      if(inputs->yEnd + remainder < header.ny) inputs->yEnd += remainder;
      else
      {
        inputs->yEnd = header.ny - 1;
      }

      ss << "New yEnd " << inputs->yEnd << std::endl;
    }

    voxelMax[1] = inputs->yEnd;
    ss << "xStart=" << inputs->xStart << " " << "xEnd=" << inputs->xEnd << std::endl;
    ss << "yStart=" << inputs->yStart << " " << "yEnd=" << inputs->yEnd << std::endl;
    /************************************************************/
  }
  else
  {
    inputs->xStart = 0;
    inputs->yStart = 0;
    inputs->zStart = 0;

    inputs->xEnd = header.nx - 1;
    inputs->yEnd = header.ny - 1;
    inputs->zEnd = header.nz - 1;
  }

  sinogram->N_r = voxelMax[0] - voxelMin[0] + 1;
  sinogram->N_t = voxelMax[1] - voxelMin[1] + 1;

  sinogram->delta_r = 1.0;
  sinogram->delta_t = 1.0;

  FEIHeader* feiHeaders = header.feiHeaders;
  if (feiHeaders != NULL)
  {
    sinogram->delta_r = feiHeaders[0].pixelsize * 1.0e9;
    sinogram->delta_t = feiHeaders[0].pixelsize * 1.0e9;
  }

  // Clear out the vector as we are going to build it up
  inputs->goodViews.resize(0,0);
  int jStart = 0;
  bool addToGoodView = false;
  for (int i = voxelMin[2]; i <= voxelMax[2]; ++i )
  {
    addToGoodView = true;
    for(size_t j = jStart; j < inputs->excludedViews.size(); ++j)
    {
      if (inputs->excludedViews[j] == i)
      {
        addToGoodView = false;
      }
    }
    if (addToGoodView == true)
    {
     // std::cout << "Adding View Index: " << i << " To goodViews vector" << std::endl;
      inputs->goodViews.push_back(i);
    }
  }

  // The number of views is the size of the vector
  sinogram->N_theta = inputs->goodViews.size();

  // Read the subvolume of the MRC file which may contain extra views
  err = reader->read(inputs->sinoFile, voxelMin, voxelMax);
  if (err < 0)
  {
    FREE_FEI_HEADERS( header.feiHeaders )
    setErrorMessage("Error Code from Reading MRC File");
    setErrorCondition(err);
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  // This data is read as a Z,Y,X array where X is the fastest moving variable and Z is the slowest
  int16_t* data = reinterpret_cast<int16_t*>(reader->getDataPointer());

  //Allocate a 3-D matrix to store the singoram in the form of a N_y X N_theta X N_x  matrix
  // Here in the actual data, Z is the slowest, then X, then Y (The Fastest) so we
  // will need to "rotate" the data in the XY plane when copying from the MRC read data into our
  // own array.
//  sinogram->counts=(DATA_TYPE***)get_3D(sinogram->N_theta,
//                                        inputs->xEnd - inputs->xStart+1,
//                                        inputs->yEnd - inputs->yStart+1,
//                                        sizeof(DATA_TYPE));
  size_t dims[3] = {sinogram->N_theta,
  inputs->xEnd - inputs->xStart+1,
  inputs->yEnd - inputs->yStart+1};
  sinogram->counts = RealVolumeType::New(dims, "Sinogram.counts");


	//If the bright field image is included initialize space for it
	/*if(inputs->BrightFieldFile != NULL)
	{
	size_t dims[3] = {sinogram->N_theta,
	inputs->xEnd - inputs->xStart+1,
	inputs->yEnd - inputs->yStart+1};
	sinogram->counts_BF = RealVolumeType::New(dims);
	sinogram->counts_BF->setName("Sinogram.counts_BrightField");
	}*/


  sinogram->angles.resize(sinogram->N_theta);
  FEIHeader* fei = NULL;
  for (uint16_t z = 0; z < sinogram->N_theta; z++)
  {
    int dataZOffset = inputs->goodViews[z] - voxelMin[2];
    // Copy the value of the tilt angle into the inputs->angles vector
    if (NULL != header.feiHeaders) {
      int offset = inputs->goodViews[z];
      fei = &(header.feiHeaders[offset]);
      if (inputs->tiltSelection == SOC::A_Tilt) {
        sinogram->angles[z] = -fei->a_tilt;
      }
      else if (inputs->tiltSelection == SOC::B_Tilt)
      {
        sinogram->angles[z] = -fei->b_tilt;
      }
    }
   // std::cout << "data_z_index: " << inputs->goodViews[z] << "  dataZOffset: " << dataZOffset << "   counts offset: " << z << std::endl;
    for (uint16_t y = 0; y < sinogram->N_t; y++)
    {
      for (uint16_t x = 0; x < sinogram->N_r; x++)
      {
        size_t index = (dataZOffset * sinogram->N_r * sinogram->N_t) + (y * sinogram->N_r) + x;
        sinogram->counts->setValue(data[index], z, x, y);
      }
    }
  }

  // Clean up all the memory associated with the MRC Reader
  reader->setDeleteMemory(true);
  reader = MRCReader::NullPointer();
  FREE_FEI_HEADERS( header.feiHeaders )


//  sinogram->N_theta = TotalNumMaskedViews;
//  sinogram->N_r = (input->xEnd - input->xStart+1);
//  sinogram->N_t = (input->yEnd - input->yStart+1);
  sinogram->R0 = -(sinogram->N_r*sinogram->delta_r)/2;
  sinogram->RMax = (sinogram->N_r*sinogram->delta_r)/2;
  sinogram->T0 =  -(sinogram->N_t*sinogram->delta_t)/2;
  sinogram->TMax = (sinogram->N_t*sinogram->delta_t)/2;


  ss << "Size of the Masked Sinogram N_r =" << sinogram->N_r << " N_t = "<< sinogram->N_t
      << " N_theta=" << sinogram->N_theta << std::endl;

  if(getVerbose())
  {

    //check sum calculation
    for (uint16_t i = 0; i < sinogram->N_theta; i++)
    {
      sum = 0;
      for (uint16_t j = 0; j < sinogram->N_r; j++)
      {
        for (uint16_t k = 0; k < sinogram->N_t; k++)
        {
          sum += sinogram->counts->getValue(i, j, k);
        }
      }
      ss << "Sinogram Checksum " << i << ":" << sum << std::endl;
    }
    std::cout << ss.str() << std::endl;
  }

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the MRC Input file", 0, UpdateProgressMessage);
}
