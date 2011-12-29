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



#include "NuisanceParamWriter.h"

#include <stdio.h>

#include <string>

#include "MXA/Utilities/MXADir.h"

#include "TomoEngine/SOC/SOCConstants.h"



#define MAKE_OUTPUT_FILE(Fp, err, outdir, filename)\
    {\
    std::string filepath(outdir);\
    filepath = filepath.append(MXADir::getSeparator()).append(filename);\
    errno = 0;\
    err = 0;\
    Fp = fopen(filepath.c_str(),"wb");\
    if (Fp == NULL) { std::cout << "Error " << errno << " Opening Output file " << filepath << std::endl; err = -1; }\
    }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NuisanceParamWriter::NuisanceParamWriter() :
m_WriteBinary(true),
m_NuisanceParams(NULL)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NuisanceParamWriter::~NuisanceParamWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NuisanceParamWriter::execute()
{
  if(NULL == getTomoInputs())
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceBinWriter::TomoInputs not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if(NULL == getSinogram())
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceBinWriter::Sinogram not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if(NULL == m_NuisanceParams)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceBinWriter:: NuisanceParams not initialized Correctly");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  std::string filepath(getTomoInputs()->outputDir);
  filepath = filepath.append(MXADir::getSeparator()).append(m_FileName);

  FILE* file = fopen(filepath.c_str(), "wb");
  if(file == 0)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceBinWriter: Error opening output file for writing");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  RealArrayType::Pointer src = RealArrayType::NullPointer();
  std::string printTitle("UNKNOWN DATA");
  switch(m_DataToWrite)
  {
    case Nuisance_I_O:
      src = m_NuisanceParams->I_0;
      printTitle = "Gains";
      break;
    case Nuisance_mu:
      src = m_NuisanceParams->mu;
      printTitle = "Offsets";
      break;
    case Nuisance_alpha:
      src = m_NuisanceParams->alpha;
      printTitle = "Variances";
      break;
    default:
      break;
  }
  if(NULL == src.get() || NULL == src->d)
  {
    setErrorCondition(-1);
    setErrorMessage("NuisanceBinWriter: The array to write was NULL");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }
  if(m_WriteBinary == true)
  {
    fwrite(src->d, sizeof(DATA_TYPE), getSinogram()->N_theta, file);
  }
  else
  {
    for (uint16_t i_theta = 0; i_theta < getSinogram()->N_theta; i_theta++)
    {
      fprintf(file, "%lf\n", src->d[i_theta]);
    }
  }
  fclose(file);

  std::cout << "************* " << printTitle << " ***************" << std::endl;
  for (uint16_t i_theta = 0; i_theta < getSinogram()->N_theta; i_theta++)
  {
    printf("%d  %lf\n", i_theta, src->d[i_theta]);
  }
  std::cout << "****************************" << std::endl;

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Writing the NuisanceParameters to a Binary File", 0, UpdateProgressMessage);
}
