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

#include "MRCWriter.h"

#include <stdio.h>
#include <iomanip>

#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/IO/MXAFileWriter64.h"


#include "TomoEngine/IO/MRCHeader.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCWriter::MRCWriter()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCWriter::~MRCWriter()
{
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCWriter::writeHeader()
{
  int err = -1;


  return err;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCWriter::write()
{
  int err = -1;

  if (m_OutputFile.empty())
  {
    return -1;
  }
  MXAFileWriter64 writer(m_OutputFile);
  bool success = writer.initWriter();
  if (false == success)
  {
    return -1;
  }
  std::cout << "Writing MRC file: " << m_OutputFile << std::endl;
  MRCHeader header;
  header.nx = m_Geometry->N_x;
  header.ny = m_Geometry->N_y;
  header.nz = m_Geometry->N_z;

  header.mode = 1;

  header.nxstart = 0;
  header.nystart = 0;
  header.nzstart = 0;
  header.mx = m_Geometry->N_x;
  header.my = m_Geometry->N_y;
  header.mz = m_Geometry->N_z;

  header.xlen = m_Geometry->N_x;
  header.ylen = m_Geometry->N_y;
  header.zlen = m_Geometry->N_z;

  header.alpha = 90.0;
  header.beta = 90.0;
  header.gamma = 90.0;

  header.mapc = 1;
  header.mapr = 2;
  header.maps = 3;

  /* ** These need to be calculated from the data ** */
  header.amin = 0.0;
  header.amax = 0.0;
  header.amean = 0.0;

  header.ispg = 0;
  header.nsymbt = 0;
  /* ** Calculate the size of the extended header */
  header.next = sizeof(FEIHeader) * m_Geometry->N_z;
  header.creatid = 0;
  ::memset(&header.extra_data, 0, 30);

  header.nint = 0;
  header.nreal = 0;
  ::memset(&header.extra_data_2, 0, 20);
  header.imodStamp = 0;
  header.imodFlags = 0;
  header.idtype = 0;
  header.lens = 0;
  header.nd1 = 0;
  header.nd2 = 0;
  header.vd1 = 0;
  header.vd2 = 0;

  ::memset(&header.tiltangles, 0, 6 * sizeof(float));
  header.xorg = 0.0f;
  header.yorg = 0.0f;
  header.zorg = 0.0f;
  header.cmap[0] = 'c';
  header.cmap[1] = 'm';
  header.cmap[2] = 'a';
  header.cmap[3] = 'p';
#if defined (CMP_WORDS_BIGENDIAN)
  header.stamp[0] = 17;
  header.stamp[1] = 17;
#else
  header.stamp[0] = 68;
  header.stamp[1] = 65;
#endif
  header.rms = 0.0f;
  header.nLabels = 3;
  for(int i = 0; i < 10; ++i)
  {
    ::memset(header.labels[i], 0, 80);
  }
  snprintf(header.labels[0], 80, "Fei Company");
  snprintf(header.labels[1], 80, "Reconstruction by EIM Tomography Engine");
  snprintf(header.labels[2], 80, "EIM TomoEngine code developer by Purdue University & BlueQuartz Software");

  // Write the header
  writer.write(reinterpret_cast<char*>(&header), 1024);

  for(uint16_t i = 0; i < m_Geometry->N_z; ++i)
  {
    FEIHeader fei;
    ::memset(&fei, 0, sizeof(FEIHeader));
    fei.pixelsize = m_Geometry->LengthX;
    writer.write(reinterpret_cast<char*>(&fei), sizeof(FEIHeader));
  }

  size_t size = sizeof(uint16_t) * m_Geometry->N_x * m_Geometry->N_y;
  uint16_t* slice = (uint16_t*)(malloc(size));
  size_t index = 0;
  DATA_TYPE*** d = m_Geometry->Object->d;
  for(uint16_t z = 0; z < m_Geometry->N_z; ++z)
  {
  //  std::cout << "Writing Z=" << z << " Layer" << std::endl;
    for(uint16_t x = 0; x < m_Geometry->N_x; ++x)
    {
      for(uint16_t y = 0; y < m_Geometry->N_y; ++y)
      {
        index = (m_Geometry->N_y * m_Geometry->N_x * 0) + (m_Geometry->N_y * x) + y;
        slice[index] = (uint16_t)(d[z][x][y]);
      }
    }
    writer.write(reinterpret_cast<char*>(slice), size);
  }
  free(slice);

  std::cout << "Done Writing MRC file." << std::endl;

  return err;
}
