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


#include "NuisanceParamReader.h"

#include <stdio.h>

#include <string>

#include "MXA/Utilities/MXADir.h"





// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NuisanceParamReader::NuisanceParamReader()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NuisanceParamReader::~NuisanceParamReader()
{

}

// -----------------------------------------------------------------------------
// Read the bin files for initialize nuisance parameters from previous scale
// -----------------------------------------------------------------------------
void NuisanceParamReader::execute()
{
  if(NULL == getSinogram())
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceParamReader::Sinogram not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

//  std::string filepath(getTomoInputs()->tempDir);
//  filepath = filepath.append(MXADir::getSeparator()).append(m_FileName);

  if (m_Data->getNDims() != 1)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceParamReader: DataArray should only have 1 dimension.");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if ( *(m_Data->getDims()) != getSinogram()->N_theta)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceParamReader: DataArray does not have the correct size. It should be SingoGram->N_theta length.");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  FILE* file = fopen(m_FileName.c_str(), "rb");
  if(file == 0)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceParamReader:  Error opening output file for writing");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  size_t nItems = fread(m_Data->getPointer(), m_Data->getTypeSize(), getSinogram()->N_theta, file);
  fclose(file);

  if (nItems != getSinogram()->N_theta)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceParamReader: Error Reading Values from file");
    notify(getErrorMessage().c_str(), 0, Observable::UpdateErrorMessage);
  }
  else
  {
    setErrorCondition(0);
    setErrorMessage("");
    notify("Done Reading the NuisanceParameters", 0, UpdateProgressMessage);
  }
}
