/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
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

//--C Includes
#include <string.h>

//-- Our own header
#include "EIMImage.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EIMImage::Pointer EIMImage::NewFromSourceMosaic(EIMImage::Pointer image, bool allocateBuffer)
{
  EIMImage::Pointer p = EIMImage::New();
  if (allocateBuffer == true)
  {
    int w, h;
    image->getImagePixelDimension(w, h);
    uint8_t* u8 = p->allocateImageBuffer(w, h, image->getManageMemory());
    if (NULL == u8)
    {
      p = EIMImage::NullPointer();
    }
      
      delete u8;
      u8 = NULL;
  }
  return p;
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EIMImage::EIMImage()
{
  m_ImagePixelDimension[0] = -1; m_ImagePixelDimension[1] = -1;
  m_ManageMemory = false;
  this->m_imageBuffer = NULL;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EIMImage::~EIMImage()
{
  if (this->m_imageBuffer != NULL
      && this->m_ManageMemory == true)
  {
    delete[] this->m_imageBuffer;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EIMImage::setImageBuffer(unsigned char* value, bool manageMemory)
{
  if (this->m_imageBuffer != NULL
      && this->m_ManageMemory == true
      && value != this->m_imageBuffer )
  {
    delete[] this->m_imageBuffer;
  }
  this->m_imageBuffer = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
uint8_t* EIMImage::allocateImageBuffer(int32_t width, int32_t height, bool manageMemory)
{
  this->deallocateImageBuffer();
  size_t total = (size_t)width * (size_t)height;
  m_imageBuffer = new unsigned char[total];
  this->m_ManageMemory = manageMemory;
  this->setImagePixelDimension(width, height);

  return m_imageBuffer;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
uint8_t* EIMImage::allocateSameSizeImage(EIMImage::Pointer image)
{
  return  allocateImageBuffer(image->getImagePixelWidth(), image->getImagePixelHeight(), true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t EIMImage::initializeImageWithSourceData(int32_t width, int32_t height, uint8_t* source)
{
  this->deallocateImageBuffer();
  size_t total = width * height;
  m_imageBuffer = new unsigned char[total];
  this->m_ManageMemory = true;
  this->setImagePixelDimension(width, height);

  uint8_t* b = static_cast<uint8_t*>(::memcpy(m_imageBuffer, source, total));
  return (b == m_imageBuffer) ? 1 : -1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t EIMImage::fillImageBuffer(uint8_t val)
{
  size_t total = m_ImagePixelDimension[0] * m_ImagePixelDimension[1];
  ::memset(m_imageBuffer, val, total);
  return (NULL != m_imageBuffer) ? 1 : -1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EIMImage::deallocateImageBuffer()
{
  if (this->m_imageBuffer != NULL
      && this->m_ManageMemory == true)
  {
    delete[] this->m_imageBuffer;
  }
  m_imageBuffer = NULL;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
unsigned char* EIMImage::getImageBuffer()
{
  return m_imageBuffer;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
uint8_t* EIMImage::getPointer(size_t index)
{
  return &(m_imageBuffer[index]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EIMImage::printSelf(std::ostream& out)
{

  out << "EIMImage Properties" << std::endl;
 // out << "  Origin:                 " << _origin[0] << ", " << _origin[1] << std::endl;
  //out << "  ImageMicronSize:        " << _micronSize[0] << " x " << _micronSize[1] << std::endl;
  out << "  ImagePixelDimension:         " << m_ImagePixelDimension[0] << " x " << m_ImagePixelDimension[1] << std::endl;
  //out << "  Scaling:                " << _scaling[0] << ", " << _scaling[1] << std::endl;
  out << "  ManageMemory:           " << m_ManageMemory << std::endl;
  out << "  ImageBuffer:            " << *m_imageBuffer << std::endl;
// _intersectedTile->printSelf(out);
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EIMImage::setImagePixelDimension(EIMImage::Pointer image)
{
  int s[2];
  image->getImagePixelDimension(s);
  m_ImagePixelDimension[0] = s[0];
  m_ImagePixelDimension[1] = s[1];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t EIMImage::getImagePixelWidth()
{
  return m_ImagePixelDimension[0];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t EIMImage::getImagePixelHeight()
{
  return m_ImagePixelDimension[1];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
size_t EIMImage::getTotalPixels()
{
  return static_cast<size_t> (m_ImagePixelDimension[0] * m_ImagePixelDimension[1]);
}

