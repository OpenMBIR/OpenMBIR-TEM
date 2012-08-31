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


#include "AvizoUniformCoordinateWriter.h"

//-- C Includes
#include <stdio.h>

//-- C++ Includes
#include <iomanip>
#include <limits>

//-- MXA Includes
#include "MXA/MXA.h"
#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/LogTime.h"
#include "MXA/Common/IO/MXAFileWriter64.h"
#include "MXA/Utilities/MXAFileInfo.h"
#include "MXA/Utilities/MXADir.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AvizoUniformCoordinateWriter::AvizoUniformCoordinateWriter() :
TomoFilter()
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
AvizoUniformCoordinateWriter::~AvizoUniformCoordinateWriter()
{
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AvizoUniformCoordinateWriter::execute()
{
   int err = write();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int AvizoUniformCoordinateWriter::write()
{
 int err = -1;
  std::stringstream ss;
 // std::cout << "Avizo Output File:\n  " << m_OutputFile << std::endl;
  if (m_OutputFile.empty())
  {
      ss.str("");
      ss << "AvizoUniformCoordinateWriter: Output File was Not Set";
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
      ss << "AvizoUniformCoordinateWriter: Error opening output file for writing. '" <<
          m_OutputFile << "'";
      setErrorCondition(-1);
      setErrorMessage(ss.str());
      notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
      return err;
  }

  std::string header = generateHeader();
  writer.writeString(header);

  err = writeData(writer);

  return err;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::string AvizoUniformCoordinateWriter::generateHeader()
{
  std::stringstream ss;
  if(m_WriteBinaryFile == true)
  {
#ifdef CMP_WORDS_BIGENDIAN
    ss << "# AmiraMesh BINARY 2.1\n";
#else
    ss << "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n";
#endif
  }
  else
  {
    ss << "# AmiraMesh 3D ASCII 2.0\n";
  }
  ss << "\n";
  ss << "# Dimensions in x-, y-, and z-direction\n";
  size_t x = m_XDims[1] - m_XDims[0];
  size_t y = m_YDims[1] - m_YDims[0];
  size_t z = m_ZDims[1] - m_ZDims[0];

  TomoInputsPtr inputs = getTomoInputs();


  ss << "define Lattice " << x << " " << y << " " << z << "\n\n";

  ss << "Parameters {\n";
  ss << "     OpenMBIRParams {\n";
  ss << "         Author \"OpenMBIR\",\n";
  ss << "         DateTime \"" << tifDateTime() << "\"\n";
  ss << "     }\n";

  ss << "     Units {\n";
  ss << "         Coordinates \"microns\"\n";
  ss << "     }\n";
  ss << "     Content \"" << x << "x" << y << "x" << z << " int, uniform coordinates\",\n";
  float origin[3] = {0.0, 0.0, 0.0};

  float res[3] = {inputs->delta_xy, inputs->delta_xy, inputs->delta_xz};

  ss << "     # Bounding Box is xmin xmax ymin ymax zmin zmax\n";
  ss << "     BoundingBox " << origin[0] << " " << origin[0] + (res[0] * x);
  ss << " " << origin[1] << " " << origin[1] + (res[1] * x);
  ss << " " << origin[2] << " " << origin[2] + (res[2] * x);
  ss << ",\n";
  ss << "     CoordType \"uniform\"\n";
  ss << "}\n\n";

  ss << "Lattice { float TomoVoxels } = @1\n\n";

  ss << "# Data section follows\n";

  return ss.str();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int AvizoUniformCoordinateWriter::writeData(MXAFileWriter64 &writer)
{
  std::string start("@1\n");
  writer.writeString(start);
  if(true == m_WriteBinaryFile)
  {
      size_t dims[2] = {m_XDims[1] - m_XDims[0], m_YDims[1] - m_YDims[0]};
      FloatImageType::Pointer sliceData = FloatImageType::New(dims, "temp slice data");
      float* slice = sliceData->getPointer(0);

      size_t index = 0;
      Real_t d = 0.0;
      size_t count = 0;

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
                  count++;
                  ++index;
              }
          }
          writer.write(reinterpret_cast<char*>(slice), sizeof(float) * dims[0] * dims[1]);
      }
  }
  else
  {
      Real_t d = 0.0;
      std::stringstream ss;
      for (int z = m_ZDims[1]-1; z >= m_ZDims[0]; z--)
      {
          for (int y = m_YDims[0]; y < m_YDims[1]; ++y)
          {
              ss.str("");
              for (int x = m_XDims[0]; x < m_XDims[1]; ++x)
              {
                  d = m_Geometry->Object->getValue(z, x, y);
                  ss << d << " ";
              }
              ss << "\n";
              writer.writeString(ss.str());
          }
      }
  }
  return 1;
}
