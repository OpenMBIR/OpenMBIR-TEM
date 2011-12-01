/*
 * MRCReader.cpp
 *
 *  Created on: Nov 22, 2011
 *      Author: mjackson
 */

#include "MRCReader.h"

#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCReader::MRCReader(bool deleteMemory) :
m_Header(NULL),
m_UInt8Data(NULL),
m_Int16Data(NULL),
m_UInt16Data(NULL),
m_FloatData(NULL),
m_DeleteMemory(deleteMemory)
{
  m_Header = new MRCHeader;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCReader::MRCReader() :
m_Header(NULL),
m_UInt8Data(NULL),
m_Int16Data(NULL),
m_UInt16Data(NULL),
m_FloatData(NULL),
m_DeleteMemory(true)
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
    //  std::cout << "File Offset: " << offset << std::endl;
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


