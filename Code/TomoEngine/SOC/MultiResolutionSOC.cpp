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



#include "MultiResolutionSOC.h"

#include <iostream>

#include "MXA/Utilities/MXADir.h"

#include "TomoEngine/SOC/SOCEngine.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionSOC::MultiResolutionSOC()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionSOC::~MultiResolutionSOC()
{
}


#define PRINT_VAR(out, inputs, var)\
	out << #var << ": " << inputs->var << std::endl;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionSOC::printInputs(TomoInputsPtr inputs, std::ostream &out)
{
	out << "------------------ TomoInputs Begin ------------------" << std::endl;
	PRINT_VAR(out, inputs, NumIter);
	PRINT_VAR(out, inputs, NumOuterIter);
	PRINT_VAR(out, inputs, SigmaX);
	PRINT_VAR(out, inputs, p);
	PRINT_VAR(out, inputs, StopThreshold);
	PRINT_VAR(out, inputs, InterpFlag);
	PRINT_VAR(out, inputs, useSubvolume);
	PRINT_VAR(out, inputs, xStart);
	PRINT_VAR(out, inputs, xEnd);
	PRINT_VAR(out, inputs, yStart);
	PRINT_VAR(out, inputs, yEnd);
	PRINT_VAR(out, inputs, zStart);
	PRINT_VAR(out, inputs, zEnd);
	PRINT_VAR(out, inputs, tiltSelection);
	PRINT_VAR(out, inputs, fileXSize);
	PRINT_VAR(out, inputs, fileYSize);
	PRINT_VAR(out, inputs, fileZSize);
	PRINT_VAR(out, inputs, LengthZ);
	PRINT_VAR(out, inputs, delta_xz);
	PRINT_VAR(out, inputs, delta_xy);
	PRINT_VAR(out, inputs, targetGain);
	PRINT_VAR(out, inputs, sinoFile);
	PRINT_VAR(out, inputs, initialReconFile);
	PRINT_VAR(out, inputs, gainsInputFile);
	PRINT_VAR(out, inputs, offsetsInputFile);
	PRINT_VAR(out, inputs, varianceInputFile);

	PRINT_VAR(out, inputs, tempDir);
	PRINT_VAR(out, inputs, reconstructedOutputFile);
	PRINT_VAR(out, inputs, gainsOutputFile);
	PRINT_VAR(out, inputs, offsetsOutputFile);
	PRINT_VAR(out, inputs, varianceOutputFile);

	out << "------------------ TomoInputs End ------------------" << std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionSOC::execute()
{
  std::cout << "MultiResolutionSOC::execute" << std::endl;
  int err = 0;
  std::cout << "-- There are " << m_TomoInputs.size() << " resolutions to reconstruct." << std::endl;


  //Run the First resolution to prime the pipeline
  SOCEngine::Pointer soc = SOCEngine::New();

  TomoInputsPtr prevInputs = m_TomoInputs.at(0);
  soc->setTomoInputs(prevInputs);
  SinogramPtr sinogram = SinogramPtr(new Sinogram);
  soc->setSinogram(sinogram);

  GeometryPtr geometry = GeometryPtr(new Geometry);
  soc->setGeometry(geometry);

  ScaleOffsetParamsPtr nuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);
  soc->setNuisanceParams(nuisanceParams);


  // The first time through we only have a target Gain so set the output files
  std::stringstream ss;
  ss << prevInputs->tempDir << MXADir::Separator << prevInputs->delta_xz << "x_Resolution";
  prevInputs->tempDir = ss.str(); // Reset the output Directory to a sub directory

  //Make sure the directory is created:
  bool success = MXADir::mkdir(prevInputs->tempDir, true);
  if (!success)
  {
    std::cout << "Could not create path: " << prevInputs->tempDir << std::endl;
  }

  ss.str("");
  ss << prevInputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalGainParametersFile;
  prevInputs->gainsOutputFile = ss.str();

  ss.str("");
  ss << prevInputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalOffsetParametersFile;
  prevInputs->offsetsOutputFile = ss.str();

  ss.str("");
  ss << prevInputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalVariancesFile;
  prevInputs->varianceOutputFile = ss.str();

  ss.str("");
  ss << prevInputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedBinFile;
  prevInputs->reconstructedOutputFile = ss.str();


  SOCEngine::InitializeSinogram(sinogram);
  SOCEngine::InitializeGeometry(geometry);
  SOCEngine::InitializeScaleOffsetParams(nuisanceParams);

  // We need to get messages to the gui or command line
  soc->addObserver(this);

  printInputs(prevInputs, std::cout);
#if 1
  soc->execute();
#else
     FILE* f = fopen(prevInputs->gainsOutputFile.c_str(), "wb");
     fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(prevInputs->offsetsOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(prevInputs->varianceOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(prevInputs->reconstructedOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
#endif
  soc = SOCEngine::NullPointer();


  for (size_t i = 1; i < m_TomoInputs.size(); ++i)
  {
    TomoInputsPtr inputs = m_TomoInputs.at(i);
    /* Get our input files from the last resolution iteration */
    inputs->gainsInputFile = prevInputs->gainsOutputFile;
    inputs->offsetsInputFile = prevInputs->offsetsOutputFile;
    inputs->varianceInputFile = prevInputs->varianceOutputFile;
    inputs->initialReconFile = prevInputs->reconstructedOutputFile;
    inputs->InterpFlag = 1;

    /* Now set the output files for this resolution */
    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << inputs->delta_xz << "x_Resolution";
    inputs->tempDir = ss.str(); // Reset the output Directory to a sub directory
    //Make sure the directory is created:
    success = MXADir::mkdir(inputs->tempDir, true);
    if (!success)
    {
      std::cout << "Could not create path: " << inputs->tempDir << std::endl;
    }
    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalGainParametersFile;
    inputs->gainsOutputFile = ss.str();

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalOffsetParametersFile;
    inputs->offsetsOutputFile = ss.str();

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalVariancesFile;
    inputs->varianceOutputFile = ss.str();

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedBinFile;
    inputs->reconstructedOutputFile = ss.str();

    soc = SOCEngine::New();
    soc->setTomoInputs(inputs);
    sinogram = SinogramPtr(new Sinogram);
    soc->setSinogram(sinogram);

    geometry = GeometryPtr(new Geometry);
    soc->setGeometry(geometry);

    nuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);
    soc->setNuisanceParams(nuisanceParams);


    SOCEngine::InitializeSinogram(sinogram);
    SOCEngine::InitializeGeometry(geometry);
    SOCEngine::InitializeScaleOffsetParams(nuisanceParams);

    // We need to get messages to the gui or command line
    soc->addObserver(this);
    printInputs(inputs, std::cout);
#if 1
    soc->execute();
#else
     f = fopen(inputs->gainsOutputFile.c_str(), "wb");
     fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(inputs->offsetsOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(inputs->varianceOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
    f = fopen(inputs->reconstructedOutputFile.c_str(), "wb");
    fprintf(f, "Testing\n");
    fclose(f); f = NULL;
#endif


    soc = SOCEngine::NullPointer();

    prevInputs = inputs;
  }

  updateProgressAndMessage("MultiResolution SOC Complete", 100);
  setErrorCondition(err);
}
