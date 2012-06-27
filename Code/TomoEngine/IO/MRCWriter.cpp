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
#include <limits>

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
void MRCWriter::execute()
{
    write();
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
  std::stringstream ss;
//  std::cout << "MRC Output File: '" << m_OutputFile << "'" << std::endl;
  if (m_OutputFile.empty())
  {
      ss.str("");
      ss << "MRCWriter: Output File was Not Set";
      setErrorCondition(-1);
      setErrorMessage(ss.str());
      notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
      return err;
  }
  MXAFileWriter64 writer(m_OutputFile);
  bool success = writer.initWriter();

  if (false == success)
  {
      ss.str("");
      ss << "MRCWriter: Error opening output file for writing. '" <<
          m_OutputFile << "'";
      setErrorCondition(-1);
      setErrorMessage(ss.str());
      notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
      return err;
  }

//  std::cout << "Writing MRC File with Geometry of: " << std::endl;
//  std::cout << "  N_z: " << m_Geometry->N_z << std::endl;
//  std::cout << "  N_x: " << m_Geometry->N_x << std::endl;
//  std::cout << "  N_y: " << m_Geometry->N_y << std::endl;


  MRCHeader header;
  header.nx = m_Geometry->N_x;
  header.ny = m_Geometry->N_y;
  header.nz = m_Geometry->N_z;

  header.mode = 2;

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
  snprintf(header.labels[2], 80, "EIM TomoEngine code developed by Purdue University & BlueQuartz Software");

  // Write the header
 // std::cout << "  Writing Header to file at " << writer.getFilePointer64() << std::endl;
  writer.write(reinterpret_cast<char*>(&header), 1024);

 // std::cout << "  Writing Extended Header to file at " << writer.getFilePointer64() << std::endl;

  for(uint16_t i = 0; i < m_Geometry->N_z; ++i)
  {
    FEIHeader fei;
    ::memset(&fei, 0, sizeof(FEIHeader));
    fei.pixelsize = static_cast<float>(m_Geometry->LengthX);
    writer.write(reinterpret_cast<char*>(&fei), sizeof(FEIHeader));
  }
  //std::cout << "  Writing Data to file at " << writer.getFilePointer64() << std::endl;

  size_t size = sizeof(float) * m_Geometry->N_x * m_Geometry->N_y;
  float* slice = (float*)(malloc(size));
  size_t index = 0;
  Real_t d = 0.0;

#if 0
  // Find the max and min value so we can scale correctly (or at least try)
  size_t nVoxels = m_Geometry->N_z * m_Geometry->N_y * m_Geometry->N_x;
  Real_t dmax = std::numeric_limits<Real_t>::min();
  Real_t dmin = std::numeric_limits<Real_t>::max();

  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(d[i] > dmax) dmax = d[i];
    if(d[i] < dmin) dmin = d[i];
  }

  Real_t k_1 = 1.0;
  if (dmax - dmin > 0.0) {
    k_1 = std::numeric_limits<uint16_t>::max() / (dmax - dmin);
  }
#endif


#if 0
  float t;
  for (int z = m_Geometry->N_z - 1; z >= 0; z--)
  {
    for (int y = 0; y < m_Geometry->N_y; ++y)
    {
      for (int x = 0; x < m_Geometry->N_x; x++)
      {
        t = static_cast<float>(m_Geometry->Object->getValue(z, x, y));
      }
    }
  }
#endif
  Real_t dmax = std::numeric_limits<Real_t>::min();
  Real_t dmin = std::numeric_limits<Real_t>::max();
 // std::cout << "Values Written to MRC File" << std::endl;
 // float t = 0.0f;
  for (int z = m_Geometry->N_z - 1; z >= 0; z--)
  {
    index = 0;
    for (int y = 0; y < m_Geometry->N_y; ++y)
    {
      for (int x = 0; x < m_Geometry->N_x; ++x)
      {
        //index = (x * m_Geometry->N_y) + y;
        d = m_Geometry->Object->getValue(z, x, y);
        slice[index] = static_cast<float>(d);
        if(d > dmax) dmax = d;
        if(d < dmin) dmin = d;
        ++index;
      }
    }
    writer.write(reinterpret_cast<char*>(slice), size);
 //   std::cout << "  Writing Slice " << z << " to file at " << writer.getFilePointer64() << std::endl;
  }


//  std::cout << "  Min float MRC Value:" << dmin << std::endl;
//  std::cout << "  Max float MRC Value:" << dmax << std::endl;
//  std::cout << "-----------------------------" << std::endl;

  free(slice);


  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Writing the MRC Output File", 0, UpdateProgressMessage);
  err = 1;
  return err;
}
