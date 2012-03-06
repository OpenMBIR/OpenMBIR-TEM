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

  SOCEngine::InitializeTomoInputs(prevInputs);
  SOCEngine::InitializeSinogram(sinogram);
  SOCEngine::InitializeGeometry(geometry);
  SOCEngine::InitializeScaleOffsetParams(nuisanceParams);

  // We need to get messages to the gui or command line
  soc->addObserver(this);

  soc->execute();
  soc = SOCEngine::NullPointer();


  for (size_t i = 1; i < m_TomoInputs.size(); ++i)
  {
    TomoInputsPtr inputs = m_TomoInputs.at(i);
    inputs->GainsFile = prevInputs->GainsFile;
    inputs->OffsetsFile = prevInputs->OffsetsFile;
    inputs->VarianceFile = prevInputs->VarianceFile;
    inputs->InitialReconFile = prevInputs->SinoFile;

    soc = SOCEngine::New();
    soc->setTomoInputs(inputs);
    sinogram = SinogramPtr(new Sinogram);
    soc->setSinogram(sinogram);

    geometry = GeometryPtr(new Geometry);
    soc->setGeometry(geometry);

    nuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);
    soc->setNuisanceParams(nuisanceParams);

    SOCEngine::InitializeTomoInputs(inputs);
    SOCEngine::InitializeSinogram(sinogram);
    SOCEngine::InitializeGeometry(geometry);
    SOCEngine::InitializeScaleOffsetParams(nuisanceParams);

    // We need to get messages to the gui or command line
    soc->addObserver(this);

    soc->execute();
    soc = SOCEngine::NullPointer();

    prevInputs = inputs;

  }

  updateProgressAndMessage("MultiResolution SOC Complete", 100);
  setErrorCondition(err);
}
