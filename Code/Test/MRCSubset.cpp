/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Singanallur Venkatakrishnan (Purdue University)
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
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Pudue
 * Univeristy, BlueQuartz Software nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
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
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdlib.h>

#include <iostream>

#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>

//-- MXA Includes
#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/IO/MXAFileWriter64.h"
#include "MXA/Utilities/MXAFileInfo.h"
#include "MXA/Utilities/MXADir.h"


#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/IO/MRCReader.h"
#include "TomoEngine/IO/MRCWriter.h"


/**
 * @brief Parses numeric values from a delimited string into a preallocated array storage.
 * The programmer MUST know in advance how many values there will be.
 * @param values The string to be parsed
 * @param format The stdio format specifier to use (%f for floats, %d for integers
 * @param output The output location to store the parsed values
 * @return Error condition
 */
template<typename T>
int parseValues(const std::string &values, const char* format, T* output)
{
  std::string::size_type pos = values.find(",", 0);
  size_t index = 0;
  int n = sscanf(values.substr(0, pos).c_str(), format, &(output[index]));
  if(n != 1)
  {
    return -1;
  }

  ++index;
  while (pos != std::string::npos && pos != values.size() - 1)
  {
    n = sscanf(values.substr(pos + 1).c_str(), format, &(output[index]));
    pos = values.find(",", pos + 1);
    ++index;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{

  TCLAP::CmdLine cmd("", ' ', TomoEngine::Version::Complete);

  TCLAP::ValueArg<std::string> inputFile("", "inputfile", "The Input MRC File", true, "", "");
  cmd.add(inputFile);

  TCLAP::ValueArg<std::string> outputFile("", "outputfile", "Output MRC FIle", true, "", "");
  cmd.add(outputFile);

  TCLAP::ValueArg<std::string> subset("", "subset", "Subset to Reconstruct in the form xmin,xmax,ymin,ymax", true, "", "");
  cmd.add(subset);


  int subsetValues[4];
  ::memset(subsetValues, 0, 6 * sizeof(int));

  try
  {
    cmd.parse(argc, argv);

    if(subset.getValue().length() != 0)
    {
      int err = parseValues(subset.getValue(), "%d", subsetValues);
      if(err < 0)
      {
        std::cout << "Error Parsing the Subvolume Dimensions. They should be entered as --subvolume 64,128,256,80,150,280" << std::endl;
        return -1;
      }

    }
  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << " error: " << e.error() << " for arg " << e.argId() << std::endl;
    std::cout << "** Unknown Arguments. Displaying help listing instead. **" << std::endl;
    return -1;
  }

  std::string filepath = inputFile.getValue();
  MRCReader::Pointer reader = MRCReader::New(true);
  MRCHeader header;
  int err = reader->readHeader(filepath.c_str(), &header);
  if(err < 0)
  {
    std::cout << "Error reading header from file '" << inputFile.getValue() << "'" << std::endl;
    return EXIT_FAILURE;
  }

  // Get the full extended header from the file
  std::vector<uint8_t> extendedHeader = reader->extendedHeader();

  reader->printHeader(&header, std::cout);


  int voxelMin[3] =
  { 0, 0, 0 };
  int voxelMax[3] =
  { header.nx - 1, header.ny - 1, 0 };

  // Get the subset of the image as a dimension
  int nx = subsetValues[1] - subsetValues[0];
  int ny = subsetValues[3] - subsetValues[2];

  // Create a new header that we can use to write out the output MRC file
  MRCHeader outHeader;
  ::memcpy(&outHeader, &header, sizeof(MRCHeader));
  outHeader.feiHeaders = header.feiHeaders; // The FEI Headers stay the same
  outHeader.nx = nx;
  outHeader.ny = ny;
  outHeader.mx = nx;
  outHeader.my = ny;
  outHeader.xlen = nx;
  outHeader.ylen = ny;



  MXAFileWriter64 writer(outputFile.getValue());
  bool success = writer.initWriter();
  std::stringstream ss;
  if (false == success)
  {
      ss.str("");
      ss << "MRCSubset: Error opening output file for writing. '" <<
          outputFile.getValue() << "'";
      return EXIT_FAILURE;
  }

  // Write the output header
  writer.write(reinterpret_cast<char*>(&outHeader), 1024);
  std::cout << "File Pointer: " << writer.getFilePointer64() << std::endl;

  writer.writeArray( &(extendedHeader.front()), header.next);

  // Write the FEI extended headers
//  for(uint16_t i = 0; i < header.nz; ++i)
//  {
//    writer.write(reinterpret_cast<char*>( &(outHeader.feiHeaders[i])), sizeof(FEIHeader));
//  }
  std::cout << "File Pointer: " << writer.getFilePointer64() << std::endl;



  std::string path = MXAFileInfo::parentPath(outputFile.getValue());
  MXADir::mkdir(path, true);

  size_t typeSize = 0;

   switch(header.mode)
   {
     case 0:
       typeSize = 1;
       break;
     case 1:
       typeSize = 2;
       break;
     case 2:
       typeSize = 4;
       break;
     case 3:
       break;
     case 4:
       break;
     case 6:
       typeSize = 2;
       break;
     case 16:
       break;
   }


  // Allocate a region of memory large enough to hold the subset region
//   uint8_t* dataPtr = reinterpret_cast<uint8_t*>(malloc(outHeader.nx * outHeader.ny * typeSize));
//   uint8_t* dest = dataPtr;


  // Read each slice
  for (int z = 0; z < header.nz; ++z)
  {
    MRCReader::Pointer mrcTiltDataReader = MRCReader::New(true);
    std::cout << "Reading Section " << z << " out of " << header.nz << std::endl;
    voxelMin[2] = z;
    voxelMax[2] = z;
    err = mrcTiltDataReader->read(filepath, voxelMin, voxelMax);

    if(err < 0)
    {
      std::cout << "Error Code from Reading: " << err << std::endl;
      return EXIT_FAILURE;
    }
    uint8_t* inputData = reinterpret_cast<uint8_t*>(mrcTiltDataReader->getDataPointer());
    uint8_t* src = NULL;
    //Copy the data, scan line by scan line from the source to the destination
    size_t outWidth = subsetValues[1] - subsetValues[0];

    for(int row = subsetValues[2]; row < subsetValues[3]; ++row)
    {
      // Calculate the destination pointer
    //  dest = dataPtr + ((row - subsetValues[2]) * outWidth * typeSize);
      // Place the pointer to the start of the row that we want
      src = inputData + typeSize*((header.nx * row) + subsetValues[0]);


      //Copy the bytes from the source into the destination
      //::memcpy(dest, src, typeSize * outWidth);

      // Write the scan line to the output file
      writer.write(reinterpret_cast<char*>(src), typeSize * outWidth);


    }


  }

 // free(dataPtr);


  std::cout << "Done Subsetting MRC File" << std::endl;
  return EXIT_SUCCESS;
}

