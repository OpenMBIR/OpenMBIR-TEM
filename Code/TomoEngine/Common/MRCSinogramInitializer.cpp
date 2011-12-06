/*
 * MRCSinogramInitializer.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "MRCSinogramInitializer.h"

#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCSinogramInitializer::MRCSinogramInitializer() :
m_Inputs(NULL),
m_Sinogram(NULL)
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

  Sinogram* sinogram = getSinogram();
  TomoInputs* input = getInputs();
  int16_t i,j,k;
  uint16_t TotalNumMaskedViews;

  DATA_TYPE sum=0;

  MRCReader::Pointer reader = MRCReader::New(true);
  MRCHeader header;
  int err = reader->readHeader(input->SinoFile, &header);
  if (err < 0)
  {
  }

  if (header.mode != 1)
  {
    std::cout << "16 bit integers are only supported. Error at line  " << __LINE__ << " in file " << __FILE__ << std::endl;
    return;
  }

  int voxelMin[3] = {0, 0, 0};
  int voxelMax[3] = {header.nx, header.ny, header.nz-1};
  if (m_Inputs->useSubvolume == true)
  {
     voxelMin[0] = input->xStart;
     voxelMin[1] = input->yStart;
     voxelMin[2] = input->zStart;

     voxelMax[0] = input->xEnd;
     voxelMax[1] = input->yEnd;
     voxelMax[2] = input->zEnd;
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

  std::vector<bool> goodViews(header.nz, 1);
  // Lay down the mask for the views that will be excluded.
  for(std::vector<uint8_t>::size_type i = 0; i < m_Inputs->ViewMask.size(); ++i)
  {
    goodViews[m_Inputs->ViewMask[i]] = 0;
  }
  int numBadViews = 0;
  for (int i = voxelMin[2]; i <= voxelMax[2]; ++i)
  {
    if(goodViews[i] == 0)
    {
      numBadViews++;
    }
    if (feiHeaders != NULL)
    {
      sinogram->angles.push_back(-feiHeaders[i].a_tilt);
    }
  }


  TotalNumMaskedViews = header.nz - numBadViews;
  sinogram->N_theta = TotalNumMaskedViews;

  err = reader->read(input->SinoFile, voxelMin, voxelMax);
  if (err < 0)
  {
  std::cout << "Error Code from Reading: " << err << std::endl;
  return ;
  }
  int16_t* data = reinterpret_cast<int16_t*>(reader->getDataPointer());



  //Allocate a 3-D matrix to store the singoram in the form of a N_y X N_theta X N_x  matrix
  sinogram->counts=(DATA_TYPE***)get_3D(TotalNumMaskedViews,
                                        input->xEnd - input->xStart+1,
                                        input->yEnd - input->yStart+1,
                                        sizeof(DATA_TYPE));

  for (k = 0; k < sinogram->N_theta; k++)
  {
    unsigned int view_count = 0;

    for (i = 0; i < sinogram->N_t; i++)
    {
      for (j = 0; j < sinogram->N_r; j++)
      {
        //std::cout<<i<<","<<j<<","<<k<<std::endl;
        //if(Sinogram->ViewMask[k] == 1)
        sinogram->counts[k][j][i] = data[k * sinogram->N_r * sinogram->N_t + i * sinogram->N_r + j];
      }
    }
    view_count++;
  }

  // Clean up all the memory associated with the MRC Reader
  reader->setDeleteMemory(true);
  reader = MRCReader::NullPointer();

  //The normalization and offset parameters for the views
  sinogram->InitialGain=(DATA_TYPE*)get_spc(TotalNumMaskedViews, sizeof(DATA_TYPE));
  sinogram->InitialOffset=(DATA_TYPE*)get_spc(TotalNumMaskedViews, sizeof(DATA_TYPE));


//  sinogram->N_theta = TotalNumMaskedViews;
//  sinogram->N_r = (input->xEnd - input->xStart+1);
//  sinogram->N_t = (input->yEnd - input->yStart+1);
  sinogram->R0 = -(sinogram->N_r*sinogram->delta_r)/2;
  sinogram->RMax = (sinogram->N_r*sinogram->delta_r)/2;
  sinogram->T0 =  -(sinogram->N_t*sinogram->delta_t)/2;
  sinogram->TMax = (sinogram->N_t*sinogram->delta_t)/2;


  printf("Size of the Masked Sinogram N_r =%d N_t = %d N_theta=%d\n",sinogram->N_r,sinogram->N_t,sinogram->N_theta);

      //check sum calculation
  for(i=0;i<sinogram->N_theta;i++)
  {
    sum=0;
    for(j=0;j<sinogram->N_r;j++)
      for(k=0;k<sinogram->N_t;k++)
       sum+=sinogram->counts[i][j][k];
    printf("Sinogram Checksum %f\n",sum);
  }
  //end ofcheck sum


  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the MRC Input file", 0, UpdateProgressMessage);

}
