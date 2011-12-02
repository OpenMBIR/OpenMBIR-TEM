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
RawGeometryWriter::RawGeometryWriter(Geom* g) :
    m_Geometry(g)
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
int RawGeometryWriter::writeFile(const std::string &filepath)
{
  int err = 0;
  FILE* Fp = fopen(filepath.c_str(), "wb");
  if(NULL == Fp)
  {
    std::cout << "Fp: " << Fp << " errno: " << errno << std::endl;
    return -1;
  }
  DATA_TYPE buffer;

	
  for (int i = 0; i < m_Geometry->N_y; ++i)
  {
	     for (int j = 0; j < m_Geometry->N_x; ++j)
		 {
		for (int k = 0; k < m_Geometry->N_z; k++)
		{
			//	std::cout << k << std::endl;
   
        buffer = m_Geometry->Object[k][j][i];
        fwrite(&buffer, sizeof(DATA_TYPE), 1, Fp);
      }
    }
  }

  fclose(Fp);
  return err;
}
