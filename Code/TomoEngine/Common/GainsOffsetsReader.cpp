/*
 * GainsOffsetsReader.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */
#include "GainsOffsetsReader.h"
#include "TomoEngine/Common/allocate.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsReader::GainsOffsetsReader()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsReader::~GainsOffsetsReader()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsReader::execute()
{
  notify("GainsOffsetsReader Starting", 0, UpdateProgressMessage);
  Sinogram* sinogram = getSinogram();
  TomoInputs* inputs = getInputs();


  //The normalization and offset parameters for the views
  sinogram->InitialGain=(DATA_TYPE*)get_spc(sinogram->N_theta, sizeof(DATA_TYPE));
  sinogram->InitialOffset=(DATA_TYPE*)get_spc(sinogram->N_theta, sizeof(DATA_TYPE));



  FILE* Fp = NULL;

  double buffer = 0;
  Fp = fopen(inputs->GainsOffsetsFile.c_str(), "r"); //This file contains the Initial unscatterd counts and background scatter for each view
  //Fp=fopen("/Users/singanallurvenkatakrishnan/Desktop/Work/Tomography/TomoSoftware/HAADFSTEM/C-Code/Data/ConvergedGainOffsetParamsOuter45_AMConstraint.bin", "r");
  if(Fp != NULL)
  {
    std::cout << "------------Initial Gains-----------" << std::endl;
    for (int i = 0; i < sinogram->N_theta; i++)
    {

      fread( &buffer, sizeof(double), 1, Fp);
      sinogram->InitialGain[i] = buffer;
      std::cout << "Tilt: " << i << "  Gain: "<< sinogram->InitialGain[i] << std::endl;
    }
    std::cout << "------------Initial Offsets-----------" << std::endl;
    for (uint16_t i = 0; i < sinogram->N_theta; i++)
    {
      fread( &buffer, sizeof(double), 1, Fp);
      sinogram->InitialOffset[i] = buffer;
      std::cout << "Tilt: " << i << "  Offset: "<< sinogram->InitialOffset[i] << std::endl;
    }
    fclose(Fp);
    setErrorCondition(0);
    setErrorMessage("");
    notify("Done Reading the Gains and Offsets Input file", 0, UpdateProgressMessage);
  }
  else
  {
    setErrorCondition(-1);
    setErrorMessage("Could not open GainsOffsets File");
    notify(getErrorMessage().c_str(), 100, UpdateErrorMessage);
    return;
  }
}
