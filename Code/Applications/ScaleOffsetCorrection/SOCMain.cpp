/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
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
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "TomoEngine/TomoEngine.h"

#include <stdlib.h>
#include <string.h>

#include <string>
#include <iostream>


// MXA Includes
#include "MXA/Utilities/MXADir.h"

#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/SOC/SOCStructures.h"
#include "TomoEngine/SOC/SOCEngine.h"
#include "SOCArgsParser.h"

int main(int argc, char **argv)
{
  std::cout << "Starting ScaleOffsetCorrection Version " << TomoEngine::Version::Complete << std::endl;

  TomoInputsPtr inputs = TomoInputsPtr(new TomoInputs);
  SOCEngine::InitializeTomoInputs(inputs);

	TomoInputsPtr bf_inputs = TomoInputsPtr(new TomoInputs);
	SOCEngine::InitializeTomoInputs(bf_inputs);

  SOCArgsParser argParser;
  int err = argParser.parseArguments(argc, argv, inputs.get(), bf_inputs.get());
  if(err < 0)
  {
    std::cout << "Error Parsing the arguments." << std::endl;
    return EXIT_FAILURE;
  }

  argParser.printInputs(inputs, std::cout);
  argParser.printInputs(bf_inputs, std::cout);

#if 0
  char path1[MAXPATHLEN]; // This is a buffer for the text
  ::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
  getcwd(path1, MAXPATHLEN);
  std::cout << "Current Working Directory: " << path1 << std::endl;
#endif
  // Make sure the output directory is created if it does not exist
  if(MXADir::exists(inputs->tempDir) == false)
  {
    std::cout << "Output Directory '" << inputs->tempDir << "' does NOT exist. Attempting to create it." << std::endl;
    if(MXADir::mkdir(inputs->tempDir, true) == false)
    {
      std::cout << "Error creating the output directory '" << inputs->tempDir << "'\n   Exiting Now." << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Output Directory Created." << std::endl;
  }
#if 0
  ::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
  getcwd(path1, MAXPATHLEN);
  std::cout << "Current Working Directory: " << path1 << std::endl;
#endif

  // Create these variables so we
  SinogramPtr sinogram = SinogramPtr(new Sinogram);
  SinogramPtr bf_sinogram = SinogramPtr(new Sinogram);
  GeometryPtr geometry =  GeometryPtr(new   Geometry);
  ScaleOffsetParamsPtr nuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);

  SOCEngine::InitializeSinogram(sinogram);
  SOCEngine::InitializeGeometry(geometry);
  SOCEngine::InitializeScaleOffsetParams(nuisanceParams);

  SOCEngine::InitializeSinogram(bf_sinogram);

  SOCEngine::Pointer engine = SOCEngine::New();
  engine->setTomoInputs(inputs);
  engine->setSinogram(sinogram);
  engine->setGeometry(geometry);
  engine->setNuisanceParams(nuisanceParams);
	engine->setBFTomoInputs(bf_inputs);
	engine->setBFSinogram(bf_sinogram);
	Observer observer;
	engine->addObserver(&observer);

  // Run the reconstruction
  engine->execute();
  if(engine->getErrorCondition() < 0)
  {
    std::cout << "Error Reconstructing the Data" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Completed SOC Run" << std::endl;

  return EXIT_SUCCESS;
}

