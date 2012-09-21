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

//-- C Includes
#include <stdio.h>

//-- C++ Includes
#include <iomanip>
#include <limits>

//-- MXA Includes
#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/IO/MXAFileWriter64.h"




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCWriter::MRCWriter() :
TomoFilter(),
m_MRCHeader(NULL)
{
m_XDims[0] = 0;
m_XDims[1] = 0;
m_YDims[0] = 0;
m_YDims[1] = 0;
m_ZDims[0] = 0;
m_ZDims[1] = 0;
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
void MRCWriter::initializeMRCHeader(MRCHeader* header)
{

  header->nx = m_Geometry->N_x;
  header->ny = m_Geometry->N_y;
  header->nz = m_Geometry->N_z;

  header->mode = 2;

  header->nxstart = 0;
  header->nystart = 0;
  header->nzstart = 0;
  header->mx = m_Geometry->N_x;
  header->my = m_Geometry->N_y;
  header->mz = m_Geometry->N_z;

  header->xlen = m_Geometry->N_x;
  header->ylen = m_Geometry->N_y;
  header->zlen = m_Geometry->N_z;

  header->alpha = 90.0;
  header->beta = 90.0;
  header->gamma = 90.0;

  header->mapc = 1;
  header->mapr = 2;
  header->maps = 3;

  /* ** These need to be calculated from the data ** */
  header->amin = 0.0;
  header->amax = 0.0;
  header->amean = 0.0;

  header->ispg = 0;
  header->nsymbt = 0;
  /* ** Calculate the size of the extended header */
  header->next = sizeof(FEIHeader) * m_Geometry->N_z;
  header->creatid = 0;
  ::memset(&header->extra_data, 0, 30);

  header->nint = 0;
  header->nreal = 0;
  ::memset(&header->extra_data_2, 0, 20);
  header->imodStamp = 0;
  header->imodFlags = 0;
  header->idtype = 0;
  header->lens = 0;
  header->nd1 = 0;
  header->nd2 = 0;
  header->vd1 = 0;
  header->vd2 = 0;

  ::memset(&header->tiltangles, 0, 6 * sizeof(float));
  header->xorg = 0.0f;
  header->yorg = 0.0f;
  header->zorg = 0.0f;
  header->cmap[0] = 'M';
  header->cmap[1] = 'A';
  header->cmap[2] = 'P';
  header->cmap[3] = ' ';
#if defined (CMP_WORDS_BIGENDIAN)
  header->stamp[0] = 17;
  header->stamp[1] = 17;
#else
  header->stamp[0] = 68;
  header->stamp[1] = 65;
#endif
  header->stamp[2] = 0;
  header->stamp[3] = 0;
  header->rms = 0.0f;
  header->nLabels = 3;
  for(int i = 0; i < 10; ++i)
  {
    ::memset(header->labels[i], 0, 80);
  }
  snprintf(header->labels[0], 80, "Fei Company (C) Copyright 2003");
  snprintf(header->labels[1], 80, "Reconstruction by OpenMBIR");
  snprintf(header->labels[2], 80, "OpenMBIR code developed by Purdue University & BlueQuartz Software");

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCWriter::write()
{
  int err = -1;
  std::stringstream ss;
 // std::cout << "MRC Output File:\n  " << m_OutputFile << std::endl;
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
  MRCHeader mrcHeader; // Put one on the stack in case the programmer supplied a NULL pointer
  ::memset(&mrcHeader, 0, 1024);
  bool detachHeaderReference = false;
  if (NULL == m_MRCHeader)
  {
    m_MRCHeader = &mrcHeader;
    initializeMRCHeader(m_MRCHeader);
    detachHeaderReference = true;
  }
  // Update the header with our Dimensions
  m_MRCHeader->nx = (m_XDims[1] - m_XDims[0]);
  m_MRCHeader->ny = (m_YDims[1] - m_YDims[0]);
  m_MRCHeader->nz = (m_ZDims[1] - m_ZDims[0]);
  m_MRCHeader->mx = m_MRCHeader->nx;
  m_MRCHeader->my = m_MRCHeader->ny;
  m_MRCHeader->mz = m_MRCHeader->nz;

  m_MRCHeader->xlen = m_MRCHeader->nx;
  m_MRCHeader->ylen = m_MRCHeader->ny;
  m_MRCHeader->zlen = m_MRCHeader->nz;
  m_MRCHeader->next = sizeof(FEIHeader) * m_MRCHeader->nz;

  // Write the header
  writer.write(reinterpret_cast<char*>(m_MRCHeader), 1024);
  for(uint16_t i = 0; i < m_Geometry->N_z; ++i)
  {
    FEIHeader fei;
    ::memset(&fei, 0, sizeof(FEIHeader));
    fei.pixelsize = static_cast<float>(m_Geometry->LengthX);
    writer.write(reinterpret_cast<char*>(&fei), sizeof(FEIHeader));
  }

  size_t dims[2] = {m_MRCHeader->nx, m_MRCHeader->ny};
  FloatImageType::Pointer sliceData = FloatImageType::New(dims, "temp slice data");
  float* slice = sliceData->getPointer(0);

  size_t index = 0;
  Real_t d = 0.0;
  size_t count = 0;

  float mean = 0.0;
  float dmax = std::numeric_limits<Real_t>::min();
  float dmin = std::numeric_limits<Real_t>::max();

  for (int z = m_ZDims[1]-1; z >= m_ZDims[0]; z--)
  {
    index = 0;
    for (int y = m_YDims[0]; y < m_YDims[1]; ++y)
    {
      for (int x = m_XDims[0]; x < m_XDims[1]; ++x)
      {
        //index = (x * m_Geometry->N_y) + y;
        d = m_Geometry->Object->getValue(z, x, y);
        slice[index] = static_cast<float>(d);
        mean += slice[index];
        count++;
        if(d > dmax) dmax = d;
        if(d < dmin) dmin = d;
        ++index;
      }
    }
    writer.write(reinterpret_cast<char*>(slice), sizeof(float) * dims[0] * dims[1]);
  }

  // Calculate the mean value
  mean = mean / count;
  // Update the values in the header of the file
  writer.setFilePointer64(76); // Set the position to the "amin" header entry
  writer.writeValue(&dmin);
  writer.writeValue(&dmax);
  writer.writeValue(&mean);


  if (detachHeaderReference == true)
  {
    m_MRCHeader = NULL;
  }

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Writing the MRC Output File", 0, UpdateProgressMessage);
  err = 1;
  return err;
}
