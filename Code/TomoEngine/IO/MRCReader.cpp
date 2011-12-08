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

#include "MRCReader.h"
#include <stdio.h>
#include <iomanip>

#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCReader::MRCReader(bool deleteMemory) :
m_DeleteMemory(deleteMemory),
m_Header(NULL),
m_UInt8Data(NULL),
m_Int16Data(NULL),
m_UInt16Data(NULL),
m_FloatData(NULL)
{
  m_Header = new MRCHeader;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCReader::MRCReader() :
m_DeleteMemory(true),
m_Header(NULL),
m_UInt8Data(NULL),
m_Int16Data(NULL),
m_UInt16Data(NULL),
m_FloatData(NULL)
{

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCReader::~MRCReader()
{
  if(m_DeleteMemory) {
    if (NULL != m_UInt8Data) { free(m_UInt8Data); }
    if (NULL != m_Int16Data) { free(m_Int16Data); }
    if (NULL != m_UInt16Data) { free(m_UInt16Data); }
    if (NULL != m_FloatData) { free(m_FloatData); }
  }
  delete m_Header;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCReader::readHeader(const std::string &filepath, MRCHeader* header)
{
  MXAFileReader64 reader(filepath);
   bool success = reader.initReader();
   if (false == success)
   {
     return -1;
   }
   header->feiHeaders = NULL;
   ::memset(header, 0, 1024); // Splat zeros across the entire structure
   success = reader.rawRead(reinterpret_cast<char*>(header), 1024);
   if (false == success)
   {
     return -2;
   }

   // Now read the extended header
   std::vector<uint8_t> extended_header(header->next, 0);
   success = reader.readArray( &(extended_header.front()), header->next);
   if (false == success)
   {
     return -3;
   }

   // If we have an FEI header then parse the extended header information
   std::string feiLabel(header->labels[0], 80);
   if (feiLabel.find_first_of("Fei Company") != std::string::npos)
   {
     // Allocate and copy in the data
     header->feiHeaders = reinterpret_cast<FEIHeader*>(malloc(sizeof(FEIHeader) * header->nz));
     ::memcpy(header->feiHeaders, &(extended_header.front()), sizeof(FEIHeader) * header->nz);
   }


   return 1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::string MRCReader::getLabelField(int index)
{
  if (m_Header == NULL)
  {
    return std::string();
  }
  if (index > 9)
  {
    return std::string();
  }

  std::string label(m_Header->labels[index], 80);
  return label;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCReader::read(const std::string &filepath, int* voxelMin, int* voxelMax)
{
  bool readSubVolume = false;
  MXAFileReader64 reader(filepath);
  bool success = reader.initReader();
  if (false == success)
  {
    return -1;
  }

  m_Header = new MRCHeader;
  ::memset(m_Header, 0, 1024); // Splat zeros across the entire structure
  success = reader.rawRead(reinterpret_cast<char*>(m_Header), 1024);
  if (false == success)
  {
    return -2;
  }

  // Now read the extended header
  std::vector<uint8_t> extended_header(m_Header->next, 0);
  success = reader.readArray( &(extended_header.front()), m_Header->next);
  if (false == success)
  {
    return -3;
  }

  // If we have an FEI header then parse the extended header information
  std::string feiLabel(m_Header->labels[0], 80);
  m_Header->feiHeaders = NULL;
  if (feiLabel.find_first_of("Fei Company") != std::string::npos)
  {
    // Allocate and copy in the data
    m_Header->feiHeaders = reinterpret_cast<FEIHeader*>(malloc(sizeof(FEIHeader) * m_Header->nz));
    ::memcpy(m_Header->feiHeaders, &(extended_header.front()), sizeof(FEIHeader) * m_Header->nz);
  }



  size_t nVoxels = m_Header->nx * m_Header->ny * m_Header->nz;
  if ( NULL != voxelMin && NULL != voxelMax)
  {
    // The user is requesting a sub-volume so calculate the number of voxels to be read
    nVoxels = (voxelMax[0] - voxelMin[0] + 1) * (voxelMax[1] - voxelMin[1] + 1) * (voxelMax[2] - voxelMin[2] + 1);
    readSubVolume = true;
  }
  size_t typeSize = 0;
  void* dataPtr = NULL;
  switch(m_Header->mode)
  {
    case 0:
      m_UInt8Data = (uint8_t*)(malloc(nVoxels * sizeof(uint8_t)));
      typeSize = 1;
      dataPtr = m_UInt8Data;
      break;
    case 1:
      m_Int16Data = (int16_t*)(malloc(nVoxels * sizeof(int16_t)));
      typeSize = 2;
      dataPtr = m_Int16Data;
      break;
    case 2:
      m_FloatData = (float*)(malloc(nVoxels * sizeof(float)));
      typeSize = 4;
      dataPtr = m_FloatData;
      break;
    case 3:
      break;
    case 4:
      break;
    case 6:
      m_UInt16Data = (uint16_t*)(malloc(nVoxels * sizeof(int16_t)));
      typeSize = 2;
      dataPtr = m_UInt16Data;
      break;
    case 16:
      break;
  }

  if (typeSize > 0 && dataPtr != NULL )
  {
    if (false == readSubVolume) {
      success = reader.rawRead(reinterpret_cast<char*>(dataPtr), typeSize * nVoxels);
      if(false == success)
      {
        return -4;
      }
    }
    else
    {
      // Read a subvolume.
      success = readPartialVolume(reader, reinterpret_cast<char*>(dataPtr), typeSize, nVoxels, voxelMin, voxelMax);
    }
  }
  else
  {
    return -5;
  }

  return 1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool MRCReader::readPartialVolume(MXAFileReader64 &reader, char* dataPtr,
                              size_t typeSize, size_t nVoxels,
                              int* voxelMin, int* voxelMax)
{
  int64_t fileStartPos = reader.getFilePointer64();
  int64_t offset = 0;
  std::streamsize numBytes = 0;
  bool success = false;
  for (int z = voxelMin[2]; z <= voxelMax[2]; ++z)
  {

    for (int y = voxelMin[1]; y <= voxelMax[1]; ++y)
    {

      offset = (m_Header->nx * m_Header->ny * z) + (m_Header->nx * y) + voxelMin[0];
      offset = offset * typeSize; // This gives number of bytes to the start of this set of data
    //  std::cout << "File Offset: " << offset )
      numBytes = (voxelMax[0] - voxelMin[0] + 1) * typeSize; // This gives number of bytes to read
      reader.setFilePointer64(fileStartPos + offset);
      success = reader.rawRead(dataPtr, numBytes);
      if (false == success)
      {
        return false;
      }
      dataPtr = dataPtr + numBytes; // Move the pointer forward the number of bytes that was read
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void* MRCReader::getDataPointer()
{
  if (NULL == m_Header)
  {
    return NULL;
  }
  switch(m_Header->mode)
  {
    case 0:
      return m_UInt8Data;
      break;
    case 1:
      return m_Int16Data;
      break;
    case 2:
      return m_FloatData;
      break;
    case 3:
      break;
    case 4:
      break;
    case 6:
      return m_UInt16Data;
      break;
    case 16:
      break;
  }
  return NULL;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCHeader* MRCReader::getHeader()
{
  return m_Header;
}

#define PRINT_VARIABLE(out, desc, var)\
  out << desc << "  " << var << std::endl;

#define FW 14

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCReader::printHeader(MRCHeader* h, std::ostream &out)
{
  out << "MRC Header ----------------------------------" << std::endl;
  PRINT_VARIABLE (out, "nx:" , h->nx )
  PRINT_VARIABLE (out, "ny:" , h->ny )
  PRINT_VARIABLE (out, "nz:" , h->nz )
  PRINT_VARIABLE (out, "mode:" , h->mode )
  PRINT_VARIABLE (out, "nxstart:" , h->nxstart )
  PRINT_VARIABLE (out, "nystart:" , h->nystart )
  PRINT_VARIABLE (out, "nzstart:" , h->nzstart )
  PRINT_VARIABLE (out, "mx:" , h->mx )
  PRINT_VARIABLE (out, "my:" , h->my )
  PRINT_VARIABLE (out, "mz:" , h->mz )
  PRINT_VARIABLE (out, "xlen:" , h->xlen )
  PRINT_VARIABLE (out, "ylen:" , h->ylen )
  PRINT_VARIABLE (out, "zlen:" , h->zlen )
  PRINT_VARIABLE (out, "alpha:" , h->alpha )
  PRINT_VARIABLE (out, "beta:" , h->beta )
  PRINT_VARIABLE (out, "gamma:" , h->gamma )
  PRINT_VARIABLE (out, "mapc:" , h->mapc )
  PRINT_VARIABLE (out, "mapr:" , h->mapr )
  PRINT_VARIABLE (out, "maps" , h->maps )
  PRINT_VARIABLE (out, "amin:" , h->amin )
  PRINT_VARIABLE (out, "amax:" , h->amax )
  PRINT_VARIABLE (out, "amean:" , h->amean )
  PRINT_VARIABLE (out, "ispg:" , h->ispg )
  PRINT_VARIABLE (out, "nsymbt:" , h->nsymbt )
  PRINT_VARIABLE (out, "next:" , h->next )
  PRINT_VARIABLE (out, "creatid:" , h->creatid )
  PRINT_VARIABLE (out, "extra_data:" , h->extra_data ) // print hex
  PRINT_VARIABLE (out, "nreal:" , h->nreal )
  PRINT_VARIABLE (out, "extra_data_2:" , h->extra_data_2 )
  PRINT_VARIABLE (out, "imodStamp:" , h->imodStamp )
  PRINT_VARIABLE (out, "imodFlags:" , h->imodFlags )
  PRINT_VARIABLE (out, "idtype:" , h->idtype )
  PRINT_VARIABLE (out, "lens:" , h->lens )
  PRINT_VARIABLE (out, "nd1:" , h->nd1 )
  PRINT_VARIABLE (out, "nd2:" , h->nd2 )
  PRINT_VARIABLE (out, "vd1:" , h->vd1 )
  PRINT_VARIABLE (out, "vd2:" , h->vd2 )
  PRINT_VARIABLE (out, "tiltangles:" , h->tiltangles[0] << h->tiltangles[1] << h->tiltangles[2] << h->tiltangles[3] << h->tiltangles[4] << h->tiltangles[5]   )
  PRINT_VARIABLE (out, "xorg:" , h->xorg )
  PRINT_VARIABLE (out, "yorg:" , h->yorg )
  PRINT_VARIABLE (out, "zorg:" , h->zorg )
  PRINT_VARIABLE (out, "cmap:" , h->cmap[0] << h->cmap[1] << h->cmap[2] << h->cmap[3]  )
  PRINT_VARIABLE (out, "stamp:" , h->stamp[0] << h->stamp[1] << h->stamp[2] << h->stamp[3]  )
  PRINT_VARIABLE (out, "rms:" , h->rms )
  PRINT_VARIABLE (out, "nLabels" , h->nLabels )

  for(int i = 0; i < h->nLabels; ++i) {
    PRINT_VARIABLE (out, "  Label: " << i , h->labels[i] )
  }


  if (h->feiHeaders != NULL)
  {
    char buf[FW + 1];
    buf[FW] = 0;
    std::vector<std::string> strings;
    strings.push_back("a_tilt");
    strings.push_back("b_tilt");
    strings.push_back("x_stage");
    strings.push_back("y_stage");
    strings.push_back("z_stage");
    strings.push_back("x_shift");
    strings.push_back("y_shift");
    strings.push_back("defocus");
    strings.push_back("exp_time");
    strings.push_back("mean_int");
    strings.push_back("tiltaxis");
    strings.push_back("pixelsize");
    strings.push_back("magnification");
    strings.push_back("voltage");
    for(size_t i = 0; i < strings.size(); ++i)
    {
      ::memset(buf, 32, FW);
      snprintf(buf, FW, "%s", strings[i].c_str());
      buf[strings[i].length()] = 32;
      out << buf << "\t";
    }
    out << std::endl;


   // std::cout << "a_tilt \t b_tilt \t x_stage \t y_stage \t z_stage \t x_shift \t y_shift \t defocus \t exp_time \t mean_int \t tiltaxis \t pixelsize \t magnification \t voltage" << std::endl;
    FEIHeader* fei = NULL;
    for (int i = 0; i < h->nz; ++i)
    {
      fei = &(h->feiHeaders[i]);
      out.setf(std::ios::left);
      out << std::setw(FW) << fei->a_tilt << "\t" << std::setw(FW) << fei->b_tilt << "\t"
          << std::setw(FW) << fei->x_stage << "\t" << std::setw(FW) << fei->y_stage << "\t" << std::setw(FW) << fei->z_stage << "\t"
          << std::setw(FW) << fei->x_shift << "\t" << std::setw(FW) << fei->y_shift << "\t"
          << std::setw(FW) << fei->defocus << "\t" << std::setw(FW) << fei->exp_time << "\t"
          << std::setw(FW) << fei->mean_int << "\t" << std::setw(FW) << fei->tiltaxis << "\t"
          << std::setw(FW) << fei->pixelsize << "\t" << std::setw(FW) << fei->magnification << "\t"
          << std::setw(FW) << fei->voltage << "\t" << std::endl;
    }
  }

  out << "---------------------------------------" << std::endl;
}
