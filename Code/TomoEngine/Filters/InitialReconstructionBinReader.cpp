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

#include "InitialReconstructionBinReader.h"

#include <sstream>

#include "MXA/Utilities/MXADir.h"

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
  std::stringstream ss;
  TomoInputsPtr input = getTomoInputs();
  GeometryPtr geometry = getGeometry();
//	uint16_t INTERP_FACTOR=2;
  //Read the Initial Reconstruction data into a 3-D matrix
  //If Interpolate flag is set then the input has only half the
  //number of voxels along each dimension as the output

  std::stringstream outPath;
  outPath << input->tempDir << MXADir::Separator << "UpsampledObject.bin";
  FILE* Fp2 = fopen(outPath.str().c_str(), "wb");
  FILE* Fp = fopen(input->initialReconFile.c_str(), "r");
  std::cout << "Reading Geom from File: " << input->initialReconFile << std::endl;
  if(NULL == Fp)
  {
    ss << "Error opening Input file '" << input->initialReconFile << "'";
    notify(ss.str(), 0, Observable::UpdateErrorMessage);
    return;
  }
  if(NULL == Fp2)
  {
    ss << "Error opening Output file '" << input->initialReconFile << "'";
    notify(ss.str(), 0, Observable::UpdateErrorMessage);
    return;
  }

  //If we need to interpolate
  std::cout << input->InterpFlag << std::endl;
  if(1 == input->InterpFlag)
  {
    std::cout << "Interpolating the initial input file by a factor of 2" << std::endl;
    uint16_t yEnd = geometry->N_y / 2;
    uint16_t xEnd = geometry->N_x / 2;
    uint16_t zEnd = geometry->N_z / 2;

    for (uint16_t y = 0; y < yEnd; y++)
    {
      for (uint16_t x = 0; x < xEnd; x++)
      {
        for (uint16_t z = 0; z < zEnd; z++)
        {
          Real_t buffer;
          fread((unsigned char*)(&buffer), sizeof(Real_t), 1, Fp);
          //Voxel replicate
          for (uint16_t p = 0; p < 2; p++)
          {
            for (uint16_t q = 0; q < 2; q++)
            {
              for (uint16_t r = 0; r < 2; r++)
              {
                geometry->Object->setValue(buffer, 2 * z + p, 2 * x + q, 2 * y + r);
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for (uint16_t i = 0; i < geometry->N_y; i++)
    {
      for (uint16_t j = 0; j < geometry->N_x; j++)
      {
        for (uint16_t k = 0; k < geometry->N_z; k++)
        {
          Real_t buffer;
          fread((unsigned char*)(&buffer), sizeof(Real_t), 1, Fp);
          geometry->Object->setValue(buffer, k, j, i);
//          geometry->Object->d[k][j][i] = buffer;
        }
      }
    }
  }
  double buffer = 0.0;
  fclose(Fp);
  notify("Done Reading Initial Reconstruction", 0, UpdateProgressMessage);
  notify("Writng Upsampled File....", 0, UpdateProgressMessage);
  for (uint16_t i = 0; i < geometry->N_y; i++)
  {
    for (uint16_t j = 0; j < geometry->N_x; j++)
    {
      for (uint16_t k = 0; k < geometry->N_z; k++)
      {
//        Real_t buffer = geometry->Object->d[k][j][i];
        buffer = geometry->Object->getValue(k, j, i);
        fwrite(&buffer, sizeof(double), 1, Fp2);
      }
    }
  }
  fclose(Fp2);
}

