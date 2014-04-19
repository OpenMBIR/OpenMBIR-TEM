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
#include <stdlib.h>

#include <string>
#include <iostream>
#include <vector>


#include "MBIRLib/MBIRLibVersion.h"
#include "MBIRLib/IOFilters/MRCReader.h"
#include "MBIRLib/IOFilters/MRCHeader.h"

#include "MXA/Utilities/MXAFileInfo.h"
#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/StringUtils.h"
//#include "MXA/Common/LogTime.h"


#if EIMTomo_TIFF_SUPPORT
#define _TIFF_DATA_TYPEDEFS_ 1
#include <tiffio.h>

int writeColorTiff(const std::string filename, std::vector<uint8_t> image, int width, int height,
                   const std::string imageDescription, int orientation)
{

  int err;
  TIFF* out;
  std::string dateTime;
  char software[1024];
  tsize_t area;

  if (image.size() == 0)
  {
    return -1;
  }
  out = TIFFOpen(filename.c_str(), "w");
  if (out == NULL)
  {
    printf("Could not open output file '%s' for writing.\n", filename.c_str());
    return -1;
  }

  err = 0;
  // set the basic values
  err = TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
  err = TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
  err = TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
  err = TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
  err = TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height); // 1 strip
  err = TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG); // single image plane

  dateTime = tifDateTime();
  err = TIFFSetField(out, TIFFTAG_DATETIME, dateTime.c_str());
  // String based tags
  if (filename.empty() == false)
  {
    err = TIFFSetField(out, TIFFTAG_DOCUMENTNAME, filename.c_str());
  }
  if (imageDescription.empty() == false)
  {
    err = TIFFSetField(out, TIFFTAG_IMAGEDESCRIPTION, imageDescription.c_str());
  }

  err = TIFFSetField(out, TIFFTAG_ORIENTATION, orientation);
  err = TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);


#if USE_LZW_COMPRESSION
  err = TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
  err = TIFFSetField(image, TIFFTAG_PREDICTOR, PREDICTOR_HORIZONTAL);
#else
  err = TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
#endif

  // Insert Resolution Units here if possible


  memset(software, 0, 1024);
  snprintf(software, 1024, "%s using libTif", TomoEngine_PACKAGE_COMPLETE);

  err = TIFFSetField(out, TIFFTAG_SOFTWARE, software);

  err = TIFFSetField(out, TIFFTAG_HOSTCOMPUTER, MXADATAMODEL_SYSTEM);

  // Write the information to the file
  area = (tsize_t)( width *  height * 3);
  err = TIFFWriteEncodedStrip(out, 0, &(image.front()), area);
  if (err != area)
  {
    err = -1;
  }
  else
  {
    err = 1;
  }

  (void)TIFFClose(out);
  return err;
}
#endif

///getColorCorrespondingToValue ////////////////////////////////////////////////
//
// Assumes you've already generated min and max -- the extrema for the data
// to which you're applying the color map. Then define the number of colorNodes
// and make sure there's a row of three float values (representing r, g, and b
// in a 0.0-1.0 range) for each node. Then call this method for with parameter
// val some float value between min and max inclusive. The corresponding rgb
// values will be returned in the reference-to-float parameters r, g, and b.
//
////////////////////////////////////////////////////////////////////////////////
void getColorCorrespondingTovalue(int16_t val,
                                  float& r, float& g, float& b,
                                  float max, float min)
{
  static const int numColorNodes = 4;
  float color[numColorNodes][3] =
  {
    {0.25f, 0.2549f, 0.7961f},    // blue
    {0.8274f, 0.8039f, 0.0941f},    // yellow
    {0.1803f, 0.6655f, 0.1490f},    // Green
    {1.0f, 0.0f, 0.0f}     // red
  };
  float range = max - min;
  for (int i = 0; i < (numColorNodes - 1); i++)
  {
    float currFloor = min + ((float)i / (numColorNodes - 1)) * range;
    float currCeil = min + ((float)(i + 1) / (numColorNodes - 1)) * range;

    if((val >= currFloor) && (val <= currCeil))
    {
      float currFraction = (val - currFloor) / (currCeil - currFloor);
      r = color[i][0] * (1.0f - currFraction) + color[i + 1][0] * currFraction;
      g = color[i][1] * (1.0f - currFraction) + color[i + 1][1] * currFraction;
      b = color[i][2] * (1.0f - currFraction) + color[i + 1][2] * currFraction;
    }
  }
}




//const std::string filepath("/Users/Shared/Data/TomographyData/TiO2Ps100kRun3/Run3TiO2PS100k_2.ali");
//const std::string filepath("/Users/Shared/Data/TomographyData/TiO2Ps100kRun3/Run3TiO2PS100k.mrc");

int main(int argc, char** argv)
{
  if(argc != 3)
  {
    std::cout
        << "This program needs the path to the mrc file and either a 0 or 1 to indicate if tiff images should be output. Output files will be placed in a folder next to that file."
        << std::endl;
    return 1;
  }
  std::string filepath(argv[1]);
  std::cout << "Testing file \n  " << filepath << std::endl;

  MRCReader::Pointer reader = MRCReader::New(true);
  MRCHeader header;
  int err = reader->readHeader(filepath, &header);
  if (err < 0)
  {
    std::cout << "Error reading header from file" << std::endl;
    return EXIT_FAILURE;
  }

  reader->printHeader(&header, std::cout);



  if ( *(argv[2]) == '0')
  {
    return EXIT_SUCCESS;
  }

#if EIMTomo_TIFF_SUPPORT


  // Write a folder full of tiff images
  int voxelMin[3] = {0, 0, 0};
  int voxelMax[3] = {header.nx - 1, header.ny - 1, 0};
//  int dims[3] = { (voxelMax[0] - voxelMin[0] + 1),
//                  (voxelMax[1] - voxelMin[1] + 1),
//                  (voxelMax[2] - voxelMin[2] + 1) };

  // Generate a Color Table
  float max = static_cast<float>(header.amax);
  float min = static_cast<float>(header.amin);
  int numColors = static_cast<int>((max - min) + 1);

  std::vector<unsigned char> colorTable(numColors * 3);
  float range = max - min;

  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int16_t val = static_cast<int16_t>( min + ((float)i / numColors) * range);
    getColorCorrespondingTovalue(val, r, g, b, max, min);
    colorTable[i * 3] = static_cast<unsigned char>(r * 255);
    colorTable[i * 3 + 1] = static_cast<unsigned char>(g * 255);
    colorTable[i * 3 + 2] = static_cast<unsigned char>(b * 255);
  }

  std::string path = MXAFileInfo::parentPath(filepath);
  path = path + MXADir::getSeparator() + MXAFileInfo::fileNameWithOutExtension(filepath) + "_Extracted_Images";
  MXADir::mkdir(path, true);

  // Read each slice
  for (int z = 0; z < header.nz; ++z)
  {
    MRCReader::Pointer r = MRCReader::New(true);
    std::cout << "Reading Section " << z << " out of " << header.nz << std::endl;
    voxelMin[2] = z;
    voxelMax[2] = z;
    err = r->read(filepath, voxelMin, voxelMax);

    if (err < 0)
    {
      std::cout << "Error Code from Reading: " << err << std::endl;
      return EXIT_FAILURE;
    }
    int16_t* data = reinterpret_cast<int16_t*>(r->getDataPointer());
    // Create an RGB Image
    std::vector<unsigned char> image(header.nx * header.ny * 3);

    int totalPixels = header.nx * header.ny;
    for (int idx = 0; idx < totalPixels; ++idx)
    {
      int colorIndex = data[idx] - static_cast<int>(header.amin);
      image[idx * 3] = colorTable[colorIndex * 3];
      image[idx * 3 + 1] = colorTable[colorIndex * 3 + 1];
      image[idx * 3 + 2] = colorTable[colorIndex * 3 + 2];
    }

    std::stringstream ss;

    std::string baseName = MXAFileInfo::fileNameWithOutExtension(filepath);
    ss << baseName;
    baseName = path + MXADir::getSeparator() + baseName.append("_").append(StringUtils::numToString(z));
    baseName = baseName.append(".tif");
    std::cout << "Writing File " << baseName << std::endl;
    ss << ": Tilt Index " << z;
    err = writeColorTiff(baseName, image, header.nx, header.ny,
                         ss.str(), ORIENTATION_TOPLEFT);
    if (err < 0)
    {
      std::cout << "Error Writing Tif file for slice " << z << std::endl;
      return EXIT_FAILURE;
    }

  }
#endif

  std::cout << "Done Reading File" << std::endl;
  return EXIT_SUCCESS;
}


