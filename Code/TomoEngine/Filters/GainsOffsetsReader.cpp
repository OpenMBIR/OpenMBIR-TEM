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
  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr inputs = getTomoInputs();


  std::vector<double> fileGains(inputs->fileZSize, 0);
  std::vector<double> fileOffsets(inputs->fileZSize, 0);

  FILE* Fp = NULL;
  size_t elementsRead = 0;
 // double buffer = 0;
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
    size_t dims[1] = {sinogram->N_theta};
    sinogram->InitialGain = RealArrayType::New(dims);
    sinogram->InitialOffset = RealArrayType::New(dims);

    // Copy just the values of the gains and offsets we need from the data read
    // from the file. The indices into the fileGains/fileOffsets array are stored
    // in the inputs->goodViews vector
    for (size_t i = 0; i < inputs->goodViews.size(); i++)
    {
//      std::cout << "Gains/Offsets Index to Copy: " << inputs->goodViews[i]
//      << " Into Index: " << i << std::endl;
      sinogram->InitialGain->d[i] = fileGains[inputs->goodViews[i]];
      sinogram->InitialOffset->d[i] = fileOffsets[inputs->goodViews[i]];
    }

#if 1
	  std::cout << "------------Initial Gains-----------" << std::endl;
	  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	  {
		  std::cout << "Tilt: " << i_theta<< "  Gain: "<< sinogram->InitialGain->d[i_theta]<< std::endl;
	  }
	  std::cout << "------------Initial Offsets-----------" << std::endl;
	  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	  {
		  std::cout << "Tilt: " << i_theta << "  Offset: "<< sinogram->InitialOffset->d[i_theta] << std::endl;
	  }
#endif
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
