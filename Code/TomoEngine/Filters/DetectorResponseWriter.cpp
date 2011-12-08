/*
 * DetectorResponseWriter.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: mjackson
 */

#include "DetectorResponseWriter.h"

#include <stdio.h>

#include <string>
#include <sstream>

#include "MXA/Utilities/MXADir.h"
#include "TomoEngine/SOC/SOCConstants.h"

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
  FILE* Fp = NULL;

  std::string filepath(getInputs()->outputDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::DetectorResponseFile);


  Fp = fopen(filepath.c_str(), "wb");
  if(Fp == NULL)
  {
    std::stringstream s;
    s << "Error Opening Detector Response Output file for writing. The output file path was \n  " << filepath;
    setErrorMessage(s.str());
    setErrorCondition(-1);
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  size_t numElements = getSinogram()->N_theta * DETECTOR_RESPONSE_BINS;
  size_t nEleWritten = fwrite(&(m_Response->d[0][0][0]), m_Response->getTypeSize(), numElements , Fp);
  fclose(Fp);
  if (numElements != nEleWritten)
  {
    setErrorCondition(-1);
    setErrorMessage("Error Writing the Detector Response File. The number of elements requested to be written was not honored");
    notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
    return;
  }

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Done Calculating the Detector Response", 0, UpdateProgressMessage);
}
