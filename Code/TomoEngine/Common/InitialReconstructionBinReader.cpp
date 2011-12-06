/*
 * InitialReconstructionBinReader.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "InitialReconstructionBinReader.h"

#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/EIMMath.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionBinReader::InitialReconstructionBinReader()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionBinReader::~InitialReconstructionBinReader()
{
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InitialReconstructionBinReader::execute()
{
  SuperClass::execute();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InitialReconstructionBinReader::initializeData()
{
  TomoInputs* input = getInputs();
  Geometry* geometry = getGeometry();

  //Read the Initial Reconstruction data into a 3-D matrix
  FILE* Fp=fopen(input->InitialReconFile.c_str(),"r");
  std::cout<<"Reading Geom"<<std::endl;
  if(NULL != Fp)
  {
    for (uint16_t i = 0; i < geometry->N_y; i++)
    {
      for (uint16_t j = 0; j < geometry->N_x; j++)
      {
        for (uint16_t k = 0; k < geometry->N_z; k++)
        {
          DATA_TYPE buffer;
          fread((unsigned char*)(&buffer), sizeof(DATA_TYPE), 1, Fp);
          geometry->Object[k][j][i] = buffer;
        }
      }
    }
    fclose(Fp);
    notify("Done Reading Initial Reconstruction", 0, UpdateProgressMessage);
  }
}

