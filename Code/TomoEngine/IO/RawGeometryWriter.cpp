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
#include "RawGeometryWriter.h"

#include <stdio.h>
#include <errno.h>

#include <iostream>

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RawGeometryWriter::RawGeometryWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RawGeometryWriter::~RawGeometryWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RawGeometryWriter::execute()
{
  notify("Starting RawGeometryWriter", 0, UpdateProgressMessage);
//  Sinogram* sinogram = getSinogram();
 // TomoInputs* input = getInputs();
  Geometry* geometry = getGeometry();

  FILE* Fp = fopen(getFilePath().c_str(), "wb");
  if(NULL == Fp)
  {
    std::cout << "Fp: " << Fp << " errno: " << errno << std::endl;
    setErrorCondition(-1);
    setErrorMessage("Could not open output file for writing");
    notify(getErrorMessage().c_str(), 100, UpdateErrorMessage);
    return;
  }
  DATA_TYPE buffer;

  for (uint16_t i = 0; i < geometry->N_y; ++i)
  {
    for (uint16_t j = 0; j < geometry->N_x; ++j)
    {
      for (uint16_t k = 0; k < geometry->N_z; k++)
      {
      //  std::cout << k << std::endl;
        buffer = geometry->Object->d[k][j][i];
        fwrite(&buffer, sizeof(DATA_TYPE), 1, Fp);
      }
    }
  }

  fclose(Fp);

  setErrorCondition(0);
  setErrorMessage("");
  std::string msg("Done with ");
  msg = msg.append(getNameOfClass());
  notify(msg.c_str(), 0, UpdateProgressMessage);

  return;
}
