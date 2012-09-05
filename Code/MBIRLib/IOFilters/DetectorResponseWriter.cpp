/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Singanallur Venkatakrishnan (Purdue University)
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

#include "DetectorResponseWriter.h"

#include <stdio.h>

#include <string>
#include <sstream>

#include "MXA/Utilities/MXADir.h"
#include "MXA/Common/IO/MXAFileWriter64.h"

#include "MBIRLib/Reconstruction/ReconstructionConstants.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectorResponseWriter::DetectorResponseWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectorResponseWriter::~DetectorResponseWriter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectorResponseWriter::execute()
{
 // FILE* Fp = NULL;

  std::string filepath(getTomoInputs()->tempDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::DetectorResponseFile);
//  Fp = fopen(filepath.c_str(), "wb");
//  if(Fp == NULL)
//  {
//    std::stringstream s;
//    s << "Error Opening Detector Response Output file for writing. The output file path was \n  " << filepath;
//    setErrorMessage(s.str());
//    setErrorCondition(-1);
//    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
//    return;
//  }



  MXAFileWriter64 writer(filepath);
  if (writer.initWriter() == false)
  {
    std::stringstream s;
    s << "Error Opening Detector Response Output file for writing. The output file path was \n  " << filepath;
    setErrorMessage(s.str());
    setErrorCondition(-1);
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  int nDims = m_Response->getNDims();
  size_t* dims = m_Response->getDims();
  size_t numElements = 1;
  for (int i = 0; i < nDims; ++i)
  {
    numElements = numElements * dims[i];
  }


  if (writer.writeArray(m_Response->d, numElements) == false)
  {
    std::stringstream ss;
    ss << "Error Writing the Detector Response File.";
    setErrorCondition(-1);
    setErrorMessage(ss.str());
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }



//  size_t nEleWritten = fwrite(m_Response->d, m_Response->getTypeSize(), numElements , Fp);
//  fclose(Fp);
//  if (numElements != nEleWritten)
//  {
//    std::stringstream ss;
//    ss << "Error Writing the Detector Response File. The number of elements requested to be written was not honored"
//        << ". The number requested was " << numElements << " and the number written was " << nEleWritten;
//    setErrorCondition(-1);
//    setErrorMessage(ss.str());
//    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
//    return;
//  }

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Writing the Detector Response", 0, UpdateProgressMessage);
}
