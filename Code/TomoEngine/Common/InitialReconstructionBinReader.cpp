/*
 * InitialReconstructionBinReader.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "InitialReconstructionBinReader.h"

#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionBinReader::InitialReconstructionBinReader() :
m_Inputs(NULL),
m_Sinogram(NULL),
m_Geometry(NULL)
{
 std::cout << "!!!!!____ InitialReconstructionBinReader IS LEAKING MEMORY !!!!!!!" << std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionBinReader::~InitialReconstructionBinReader()
{
}

// -----------------------------------------------------------------------------
//Finds the maximum of absolute value elements in an array
// -----------------------------------------------------------------------------
DATA_TYPE InitialReconstructionBinReader::absMaxArray(std::vector<DATA_TYPE> &Array)
{
  uint16_t i;
  DATA_TYPE max;
  max = fabs(Array[0]);
  for(i =1; i < Array.size();i++)
    if(fabs(Array[i]) > max)
      max=fabs(Array[i]);
  return max;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InitialReconstructionBinReader::execute()
{

  Sinogram* sinogram = getSinogram();
  TomoInputs* input = getInputs();
  Geometry* geometry = getGeometry();

  FILE* Fp;
  uint16_t i,j,k;
  DATA_TYPE sum=0,max;

  //Find the maximum absolute tilt angle
  max= absMaxArray(sinogram->angles);

#ifndef FORWARD_PROJECT_MODE
  input->LengthZ *= Z_STRETCH;

#ifdef EXTEND_OBJECT
  geometry->LengthX = ((sinogram->N_r * sinogram->delta_r)/cos(max*M_PI/180)) + input->LengthZ*tan(max*M_PI/180) ;
#else
  Geometry->LengthX = ((sinogram->N_r * sinogram->delta_r));
#endif //Extend object endif

#else
  Geometry->LengthX = ((Sinogram->N_r * Sinogram->delta_r));
#endif//Forward projector mode end if

//  Geometry->LengthY = (Geometry->EndSlice- Geometry->StartSlice)*Geometry->delta_xy;
  geometry->LengthY = (input->yEnd-input->yStart + 1)*sinogram->delta_t;

  geometry->N_x = ceil(geometry->LengthX/input->delta_xz);//Number of voxels in x direction
  geometry->N_z = ceil(input->LengthZ/input->delta_xz);//Number of voxels in z direction
  geometry->N_y = floor(geometry->LengthY/input->delta_xy);//Number of measurements in y direction

  printf("Geometry->Nz=%d\n",geometry->N_z);
  printf("Geometry->Nx=%d\n",geometry->N_x);
  printf("Geometry->Ny=%d\n",geometry->N_y);

  geometry->Object = (DATA_TYPE ***)get_3D(geometry->N_z, geometry->N_x, geometry->N_y, sizeof(DATA_TYPE));//Allocate space for the 3-D object
//Coordinates of the left corner of the x-z object
  geometry->x0 = -geometry->LengthX/2;
  geometry->z0 = -input->LengthZ/2;
 // Geometry->y0 = -(sinogram->N_t * sinogram->delta_t)/2 + Geometry->StartSlice*Geometry->delta_xy;
  geometry->y0 = -(geometry->LengthY)/2 ;

  //Read the Initial Reconstruction data into a 3-D matrix
  Fp=fopen(input->InitialReconFile.c_str(),"r");
/*  for(i=0;i<Geometry->N_y;i++)
    for(j=0;j<Geometry->N_x;j++)
      for(k=0;k<Geometry->N_z;k++)
  {
    fread (buffer,sizeof(DATA_TYPE),1,Fp);
    Geometry->Object[i][k][j]=*buffer;
//  printf("%f\n",Geometry->Object[i][j][k]);
  }*/

  for (i = 0; i < geometry->N_y; i++)
  {
    for (j = 0; j < geometry->N_x; j++)
    {
      for (k = 0; k < geometry->N_z; k++)
      {
        if(Fp == NULL)//If no input file has been specified or if the file does not exist just set the default values to be zero
        {
          geometry->Object[k][j][i] = 0;
        }
        else//If the iput file exists read the values
        {
          DATA_TYPE buffer;
          fread((unsigned char*)(&buffer), sizeof(DATA_TYPE), 1, Fp);
          geometry->Object[k][j][i] = buffer;
        }
      }
    }
  }

      //Doing a check sum to verify with matlab

  for(i=0;i<geometry->N_y;i++)
  {
    sum=0;
    for(j=0;j<geometry->N_x;j++)
      for(k=0;k<geometry->N_z;k++)
    {
      sum+=geometry->Object[k][j][i];
    }
    printf("Geometry check sum %f\n",sum);
  }
      //End of check sum

  fclose(Fp);





  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the Initial Reconstruction file", 0, UpdateProgressMessage);

}
