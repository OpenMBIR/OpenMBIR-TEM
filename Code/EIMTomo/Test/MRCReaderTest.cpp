/*
 * MRCReaderTest.cpp
 *
 *  Created on: Nov 22, 2011
 *      Author: mjackson
 */

#include <stdlib.h>

#include <string>
#include <iostream>
#include <vector>


#include "EIMTomo/EIMTomoVersion.h"
#include "EIMTomo/IO/MRCReader.h"
#include "EIMTomo/IO/MRCHeader.h"

#include "MXA/Utilities/MXAFileInfo.h"
#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/StringUtils.h"
#include "MXA/Common/LogTime.h"

#define _TIFF_DATA_TYPEDEFS_ 1
#include <tiffio.h>

int writeColorTiff(const std::string filename, std::vector<uint8_t> image, int width, int height,
                   const std::string imageDescription, int orientation)
{

  int err;
   TIFF *out;
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
   snprintf(software, 1024, "%s using libTif", EIMTomo_PACKAGE_COMPLETE);

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
                                   float &r, float &g, float &b,
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
      r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
      g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
      b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
    }
  }
}




const std::string filepath("/Users/Shared/Data/TomographyData/TiO2Ps100kRun3/Run3TiO2PS100k_2.ali");
//const std::string filepath("/Users/Shared/Data/TomographyData/TiO2Ps100kRun3/Run3TiO2PS100k.mrc");

int main(int argc, char **argv)
{
  std::cout << "Testing file \n  " << filepath << std::endl;

  MRCReader reader(true); // We are going to manage the data ourselves
  MRCHeader header;
  int err = reader.readHeader(filepath, &header);
  if (err < 0)
  {
    std::cout << "Error reading header from file" << std::endl;
    return EXIT_FAILURE;
  }

  reader.printHeader(&header, std::cout);

#if 0
  std::cout << "Num Cols:  " << header.nx << std::endl;
  std::cout << "Num Rows:  " << header.ny << std::endl;
  std::cout << "Num Sects: " << header.nz << std::endl;
  std::cout << "Pixel Type: " << header.mode << std::endl;

  switch(header.mode)
    {
      case 0:
        std::cout << "  Pixel Type:  Char" << std::endl;
        break;
      case 1:
        std::cout << "  Pixel Type: Signed Short" << std::endl;
        break;
      case 2:
        std::cout << "  Pixel Type: Float" << std::endl;
        break;
      case 3:
        break;
      case 4:
        break;
      case 6:
        break;
      case 16:
        break;
      default:
        std::cout << "  Unknown Pixel Type" << std::endl;
        break;
    }
  std::cout << "Min Pixel Value: " << header.amin << std::endl;
  std::cout << "Max Pixel Value: " << header.amax << std::endl;

  std::cout << "a_tilt \t b_tilt \t x_stage \t y_stage \t z_stage \t x_shift \t y_shift \t defocus \t exp_time \t mean_int \t tiltaxis \t pixelsize \t magnification \t voltage" << std::endl;
  FEIHeader* fei = NULL;
  for (int i = 0; i < header.nz; ++i)
  {
    fei = &(header.feiHeaders[i]);
    std::cout << fei->a_tilt << "'\t'" << fei->b_tilt << "'\t'"
        << fei->x_stage << "'\t'" << fei->y_stage << "'\t'" << fei->z_stage << "'\t'"
        << fei->x_shift << "'\t'" << fei->y_shift << "'\t'"
        << fei->defocus << "'\t'" << fei->exp_time << "'\t'"
        << fei->mean_int << "'\t'" << fei->tiltaxis << "'\t'"
        << fei->pixelsize << "'\t'" << fei->magnification << "'\t'"
        << fei->voltage << "'\t'" << std::endl;
  }
#endif


  int voxelMin[3] = {0,0,0};
  int voxelMax[3] = {header.nx-1, header.ny-1, 0};
  int dims[3] = { (voxelMax[0] - voxelMin[0] + 1),
                  (voxelMax[1] - voxelMin[1] + 1),
                  (voxelMax[2] - voxelMin[2] + 1) };

  // Generate a Color Table
  int max = header.amax;
  int min = header.amin;
  int numColors = (max-min) + 1;

  std::vector<unsigned char> colorTable(numColors * 3);
  int range = max - min;

  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int val = min + ((float)i / numColors) * range;
    getColorCorrespondingTovalue(val, r, g, b, max, min);
    colorTable[i*3] = r*255;
    colorTable[i*3+1] = g*255;
    colorTable[i*3+2] = b*255;
  }

  std::string path = MXAFileInfo::parentPath(filepath);
  path = path + MXADir::getSeparator() + MXAFileInfo::fileNameWithOutExtension(filepath) + "_Extracted_Images";
  MXADir::mkdir(path, true);

  // Read each slice
  for (int z = 0; z < header.nz; ++z)
  {
    MRCReader r(true);
    std::cout << "Reading Section " << z << " out of " << header.nz << std::endl;
    voxelMin[2] = z;
    voxelMax[2] = z;
    err = r.read(filepath, voxelMin, voxelMax);

    if (err < 0)
    {
      std::cout << "Error Code from Reading: " << err << std::endl;
      return EXIT_FAILURE;
    }
    int16_t* data = reinterpret_cast<int16_t*>(r.getDataPointer());
    // Create an RGB Image
    std::vector<unsigned char> image(header.nx * header.ny * 3);

    int totalPixels = header.nx * header.ny;
    for (int idx = 0; idx < totalPixels; ++idx)
    {
      int colorIndex = data[idx] - header.amin;
      image[idx * 3] = colorTable[colorIndex * 3];
      image[idx * 3 + 1] = colorTable[colorIndex * 3 + 1];
      image[idx * 3 + 2] = colorTable[colorIndex * 3 + 2];
    }


    std::string baseName = MXAFileInfo::fileNameWithOutExtension(filepath);
    baseName = path + MXADir::getSeparator() + baseName.append("_").append(StringUtils::numToString(z));
    baseName = baseName.append(".tif");
    std::cout << "Writing File " << baseName << std::endl;
    err = writeColorTiff(baseName, image, header.nx, header.ny,
                       "Tomography Slice", ORIENTATION_TOPLEFT);
    if (err < 0)
    {
      std::cout << "Error Writing Tif file for slice " << z << std::endl;
      return EXIT_FAILURE;
    }

  }






//  for(int z = voxelMin[2]; z <= voxelMax[2]; ++z)
//  {
//    for(int y = voxelMin[1]; y <= voxelMax[1]; ++y)
//    {
//      for(int x = voxelMin[0]; x <= voxelMax[0]; ++x)
//      {
//        size_t index = (dims[0] * dims[1] * z) + (dims[0] * y) + x;
//        std::cout << data[index] << " ";
//      }
//      std::cout << std::endl;
//
//    }
//  }

 // free(data);

  std::cout << "Done Reading File" << std::endl;
  return EXIT_SUCCESS;
}


