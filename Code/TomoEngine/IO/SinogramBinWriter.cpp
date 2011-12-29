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


#include "SinogramBinWriter.h"

#include <stdio.h>

#include <string>

#include "MXA/Utilities/MXADir.h"

#include "TomoEngine/SOC/SOCConstants.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SinogramBinWriter::SinogramBinWriter() :
TomoFilter()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SinogramBinWriter::~SinogramBinWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SinogramBinWriter::execute()
{
  if(NULL == getTomoInputs())
  {
    setErrorCondition(-1);
    setErrorMessage("SinogramBinWriter::TomoInputs not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if(NULL == getSinogram())
  {
    setErrorCondition(-1);
    setErrorMessage("SinogramBinWriter::Sinogram not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if(NULL == m_NuisanceParams)
  {
    setErrorCondition(-1);
    setErrorMessage("SinogramBinWriter:: NuisanceParams not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  notify("Writing Sinogram", 0, UpdateProgressMessage);

  FILE* file = NULL;

  std::string filepath(getTomoInputs()->outputDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::ReconstructedSinogramFile);\
  file = fopen(filepath.c_str(), "wb");
  if(file == 0)
  {
    setErrorCondition(-1);
    setErrorMessage("SinogramBinWriter: Error opening output file for writing");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  //Writing the final sinogram
  uint16_t thetaSize = getSinogram()->N_theta;
  uint16_t rSize = getSinogram()->N_r;
  uint16_t tSize = getSinogram()->N_t;
  DATA_TYPE value;
  for (uint16_t i_theta = 0; i_theta < thetaSize; i_theta++) // Depth
  {
    for (uint16_t i_r = 0; i_r < rSize; i_r++) // Width
    {
      for (uint16_t i_t = 0; i_t < tSize; i_t++) // Height
      {
        value = m_Data->d[i_theta][i_r][i_t];
        value *= m_NuisanceParams->I_0->d[i_theta];
        value += m_NuisanceParams->mu->d[i_theta];
        fwrite(&value, sizeof(DATA_TYPE), 1, file);
      }
    }
  }
  fclose(file);

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Writing the Sinogram to a Binary File", 0, UpdateProgressMessage);

}
