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
#ifndef _AIM_IMAGE_H_
#define _AIM_IMAGE_H_

#include <iostream>

#include "MXA/MXATypes.h"
#include "MXA/Common/MXASetGetMacros.h"


/**
* @class EIMImage EIMImage.h AIM/Common/EIMImage.h
* @brief This class represents a 2D image of type unsigned 8 bit characters.
* @author Michael A. Jackson for BlueQuartz Software
* @date Nov 5, 2009
* @version 1.0
*/
class EIMImage
{
  public:

    MXA_SHARED_POINTERS(EIMImage);
    MXA_STATIC_NEW_MACRO(EIMImage);
    MXA_TYPE_MACRO(EIMImage);
    virtual ~EIMImage();

    /**
    * @brief Creates a new EIMImage by copying all the settings from a source EIMImage
    * object and allocate the buffer at the same time. If the 'allocateBuffer' argument is
    * false then the image buffer will NOT be allocated. Note that the actual data from
    * the source mosaic is NOT copied into the new object.
    * @param image The source EIMImage to copy settings from.
    * @param allocateBuffer If false, the memory to hold the image is NOT allocated. Default=TRUE
    * @return EIMImage::Pointer object
    */
    static EIMImage::Pointer NewFromSourceMosaic(EIMImage::Pointer image, bool allocateBuffer = true);

    MXA_INSTANCE_VEC2_PROPERTY(int, ImagePixelDimension);
    const int32_t* getImagePixelDimension()
    {
      return m_ImagePixelDimension;
    }

    void setImagePixelDimension(EIMImage::Pointer image);
    int32_t getImagePixelWidth();
    int32_t getImagePixelHeight();


    void setImageBuffer(unsigned char* value, bool manageMemory = false);
    uint8_t* allocateImageBuffer(int32_t width, int32_t height, bool manageMemory = false);
    uint8_t* allocateSameSizeImage(EIMImage::Pointer image);
    unsigned char* getImageBuffer();

    uint8_t* getPointer(size_t index = 0);

    void deallocateImageBuffer();
    MXA_INSTANCE_PROPERTY(bool, ManageMemory);

    int32_t initializeImageWithSourceData(int32_t width, int32_t height, uint8_t* source);

    int32_t fillImageBuffer(uint8_t val);

    void printSelf(std::ostream& out);

    size_t getTotalPixels();

  protected:
    EIMImage();

  private:
    uint8_t*      m_imageBuffer;
    EIMImage(const EIMImage&);    // Copy Constructor Not Implemented
    void operator=(const EIMImage&);  // Operator '=' Not Implemented
};

#endif /* _AIM_IMAGE_H_ */
