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


  std::vector<double> fileGains(inputs->fileZSize, 0);
  std::vector<double> fileOffsets(inputs->fileZSize, 0);

  FILE* Fp = NULL;
  size_t elementsRead = 0;
  double buffer = 0;
  Fp = fopen(inputs->GainsOffsetsFile.c_str(), "r"); //This file contains the Initial unscatterd counts and background scatter for each view
  //Fp=fopen("/Users/singanallurvenkatakrishnan/Desktop/Work/Tomography/TomoSoftware/HAADFSTEM/C-Code/Data/ConvergedGainOffsetParamsOuter45_AMConstraint.bin", "r");
  if(Fp != NULL)
  {
    // Read all the gains in a single shot
    elementsRead = fread( &(fileGains.front()), sizeof(double), inputs->fileZSize, Fp);
    if (elementsRead != inputs->fileZSize)
    {
      setErrorCondition(-1);
      setErrorMessage("Error Reading Gains from File");
      notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
      fclose(Fp);
      return;
    }
    elementsRead = fread( &(fileOffsets.front()), sizeof(double), inputs->fileZSize, Fp);
    if (elementsRead != inputs->fileZSize)
    {
      setErrorCondition(-1);
      setErrorMessage("Error Reading Offsets from File");
      notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
      fclose(Fp);
      return;
    }

    fclose(Fp);
    // Allocate the proper amount of memory for the gains and offsets
    //The normalization and offset parameters for the views
    sinogram->InitialGain=(DATA_TYPE*)get_spc(sinogram->N_theta, sizeof(DATA_TYPE));
    sinogram->InitialOffset=(DATA_TYPE*)get_spc(sinogram->N_theta, sizeof(DATA_TYPE));

    // Copy just the values of the gains and offsets we need from the data read
    // from the file. The indices into the fileGains/fileOffsets array are stored
    // in the inputs->goodViews vector
    for (size_t i = 0; i < inputs->goodViews.size(); i++)
    {
//      std::cout << "Gains/Offsets Index to Copy: " << inputs->goodViews[i]
//      << " Into Index: " << i << std::endl;
      sinogram->InitialGain[i] = fileGains[inputs->goodViews[i]];
      sinogram->InitialOffset[i] = fileOffsets[inputs->goodViews[i]];
    }

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
