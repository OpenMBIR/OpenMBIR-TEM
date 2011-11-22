/*
 * MRCReader.cpp
 *
 *  Created on: Nov 22, 2011
 *      Author: mjackson
 */

#include "MRCReader.h"

#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/IO/MXAFileReader64.h"


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
  m_Header = new MRCHeader;
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
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int MRCReader::read(const std::string &filepath)
{

  MXAFileReader64 reader(filepath);
  bool success = reader.initReader();
  if (false == success)
  {
    return -1;
  }

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

  size_t nVoxels = m_Header->nx * m_Header->ny * m_Header->nz;
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

  if (typeSize > 0 && dataPtr != NULL) {
    success = reader.rawRead(reinterpret_cast<char*>(dataPtr), typeSize * nVoxels);
    if(false == success)
    {
      return -4;
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
void* MRCReader::getDataPointer()
{
  void* dataPtr = NULL;
  size_t typeSize = 0;
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
