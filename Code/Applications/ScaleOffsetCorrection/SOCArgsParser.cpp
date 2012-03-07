/* ============================================================================
 * Copyright (c) 2011, Singanallur Venkatakrishnan <svenkata@purdue.edu>
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
 * Neither the name of Singanallur Venkatakrishnan , Purdue University nor the
 * names of its contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
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
#include "SOCArgsParser.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>

#include "MXA/Utilities/MXADir.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/Common/EIMTime.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/IO/MRCReader.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/SOC/SOCConstants.h"

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

/**
 * @brief Parses unknown number of numeric values from a delimited string and places
 * the values into the output variable.
 * @param values The string to be parsed
 * @param format The stdio format specifier to use (%f for floats, %d for integers
 * @param output The output location to store the parsed values
 * @return Error condition
 */
template<typename T>
int parseUnknownArray(const std::string &values, const char* format, std::vector<T> &output)
{
  std::string::size_type pos = values.find(",", 0);
  size_t index = 0;
  T t;
  int n = sscanf(values.substr(0, pos).c_str(), format, &t);
  if(n != 1)
  {
    return -1;
  }
  output.push_back(t);

  ++index;
  while (pos != std::string::npos && pos != values.size() - 1)
  {
    n = sscanf(values.substr(pos + 1).c_str(), format, &(t));
    output.push_back(t);
    pos = values.find(",", pos + 1);
    ++index;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SOCArgsParser::SOCArgsParser()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SOCArgsParser::~SOCArgsParser()
{

}

int SOCArgsParser::parseArguments(int argc, char **argv, TomoInputs* inputs, TomoInputs* bf_inputs)
{
  if(NULL == inputs)
  {
    printf("The CommandLineInputspointer was null. Returning early.\n");
    return -1;
  }

  TCLAP::CmdLine cmd("", ' ', TomoEngine::Version::Complete);

  TCLAP::ValueArg<std::string> in_sinoFile("s", "sinofile", "The Sinogram File", true, "", "");
  cmd.add(in_sinoFile);

  TCLAP::ValueArg<std::string> in_InitialRecon("i", "ini_recon", "Initial Reconstruction to initialize algorithm", false, "", "");
  cmd.add(in_InitialRecon);

  //Normalizing Bright Field Scan

  TCLAP::ValueArg<std::string> in_BrightField("", "brightfield", "Acquired Bright Field tilt series", false, "", "");
  cmd.add(in_BrightField);

  // TCLAP::ValueArg<std::string> in_GainsOffsets("", "gains_offsets", "Initial Gains and Offsets to use.", false, "", "");
  // cmd.add(in_GainsOffsets);

  TCLAP::ValueArg<double> in_TargetGain("", "TargetGain", "Target Gain for unscattered electrons", true, 1, "");
  cmd.add(in_TargetGain);

  TCLAP::ValueArg<std::string> in_Gains("", "gains", "Initial Gains to use.", false, "", "");
  cmd.add(in_Gains);

  TCLAP::ValueArg<std::string> in_Offsets("", "offsets", "Initial Offsets to use.", false, "", "");
  cmd.add(in_Offsets);

  TCLAP::ValueArg<double> in_DefaultOffset("", "defaultOffset", "Default offset for all Tilts", false, 0.0, "");
  cmd.add(in_DefaultOffset);

  TCLAP::ValueArg<std::string> in_Variance("", "variance", "Initial Variance to use.", false, "", "");
  cmd.add(in_Variance);

  //Whether to interpolate initial file or not
  TCLAP::ValueArg<uint8_t> in_InterpFlag("", "interpolate", "Iterpolate Initial Reconstruction", false, 0, "");
  cmd.add(in_InterpFlag);

  TCLAP::ValueArg<std::string> in_outputFile("o", "outputfile", "The Output File", true, "", "");
  cmd.add(in_outputFile);

  TCLAP::ValueArg<std::string> in_outputDir("", "outdir", "The Output dir", true, ".", ".");
  cmd.add(in_outputDir);

  TCLAP::ValueArg<double> in_stopThreshold("T", "stopThreshold", "Stopping Threshold for inner loop", false, .009, "");
  cmd.add(in_stopThreshold);

  TCLAP::ValueArg<std::string> subvolume("", "subvolume", "SubVolume to Reconstruct in the form xmin,ymin,zmin,xmax,ymax,zmax", false, "", "");
  cmd.add(subvolume);

  TCLAP::ValueArg<int> in_numIter("n", "inner", "Number of Inner Iterations", true, 0, "0");
  cmd.add(in_numIter);
  TCLAP::ValueArg<double> in_sigmaX("l", "sigmax", "Sigma X Value", true, 1.0, "1.0");
  cmd.add(in_sigmaX);

  TCLAP::ValueArg<double> in_markov("m", "mrf", "Markov Random Field Parameter", true, 0.0, "0.0");
  cmd.add(in_markov);

  //TCLAP::ValueArg<std::string> InitialParameters("g", "initp", "InitialParameters", false, "", "");
  // cmd.add(InitialParameters);

  TCLAP::ValueArg<int> NumOuterIter("O", "outer", "NumOuterIter", true, 1, "1");
  cmd.add(NumOuterIter);

  TCLAP::ValueArg<DATA_TYPE> xz_size("", "xz_size", "Size in nm of output pixel xz plane", true, 1, "1");
  cmd.add(xz_size);
  TCLAP::ValueArg<DATA_TYPE> xy_size("", "xy_size", "Size in nm of output pixel xy plane", true, 1, "1");
  cmd.add(xy_size);

  TCLAP::ValueArg<DATA_TYPE> thickness("", "thickness", "Thickness of sample in nm", true, 0, "");
  cmd.add(thickness);

  TCLAP::ValueArg<std::string> viewMask("", "exclude_views", "Comma separated list of tilts to exclude by index", false, "", "");
  cmd.add(viewMask);

  TCLAP::ValueArg<unsigned int> tiltSelection("", "tilt_selection", "Which Tilt Values to use from file. Default is 'A'", true, SOC::A_Tilt, "0");
  cmd.add(tiltSelection);

  if(argc < 2)
  {
    std::cout << "Scale Offset Correction Command Line Version " << cmd.getVersion() << std::endl;
    std::vector<std::string> args;
    args.push_back(argv[0]);
    args.push_back("-h");
    cmd.parse(args);
    return -1;
  }

  try
  {
    cmd.parse(argc, argv);
    inputs->reconstructedOutputFile = in_outputDir.getValue() + MXADir::getSeparator() + in_outputFile.getValue();
    inputs->sinoFile = in_sinoFile.getValue();
    inputs->initialReconFile = in_InitialRecon.getValue();
    inputs->InterpFlag = in_InterpFlag.getValue();

    bf_inputs->sinoFile = in_BrightField.getValue();


    inputs->gainsInputFile = in_Gains.getValue();
    inputs->offsetsInputFile = in_Offsets.getValue(); //Offset
    inputs->varianceInputFile = in_Variance.getValue(); //variance

    inputs->gainsOutputFile = inputs->tempDir + ScaleOffsetCorrection::FinalGainParametersFile;
    inputs->offsetsOutputFile = inputs->tempDir + ScaleOffsetCorrection::FinalOffsetParametersFile;
    inputs->varianceOutputFile = inputs->tempDir + ScaleOffsetCorrection::FinalVariancesFile;

    inputs->tempDir = in_outputDir.getValue();
    inputs->NumIter = in_numIter.getValue();
    inputs->SigmaX = in_sigmaX.getValue();
    inputs->p = in_markov.getValue();
    inputs->NumOuterIter = NumOuterIter.getValue();
    inputs->StopThreshold = in_stopThreshold.getValue();
    inputs->targetGain = in_TargetGain.getValue();

    int subvolumeValues[6];
    ::memset(subvolumeValues, 0, 6 * sizeof(int));
    inputs->useSubvolume = false;

    if(subvolume.getValue().length() != 0)
    {
      int err = parseValues(subvolume.getValue(), "%d", subvolumeValues);
      if(err < 0)
      {
        std::cout << "Error Parsing the Subvolume Dimensions. They should be entered as --subvolume 64,128,256,80,150,280" << std::endl;
        return -1;
      }
      inputs->useSubvolume = true;
      inputs->xStart = subvolumeValues[0];
      inputs->xEnd = subvolumeValues[3];
      inputs->yStart = subvolumeValues[1];
      inputs->yEnd = subvolumeValues[4];
      inputs->zStart = subvolumeValues[2];
      inputs->zEnd = subvolumeValues[5];

      //Do same for BF sinogram as well

      bf_inputs->useSubvolume = true;
      bf_inputs->xStart = subvolumeValues[0];
      bf_inputs->xEnd = subvolumeValues[3];
      bf_inputs->yStart = subvolumeValues[1];
      bf_inputs->yEnd = subvolumeValues[4];
      bf_inputs->zStart = subvolumeValues[2];
      bf_inputs->zEnd = subvolumeValues[5];
    }

    for(int i = 0; i < argc; ++i)
    {
      if ( ::strncmp(argv[i], "--defaultOffset", 15) == 0)
      {
        inputs->useDefaultOffset = true;
        inputs->defaultOffset = in_DefaultOffset.getValue();
      }
    }

    inputs->delta_xz = xz_size.getValue();
    inputs->delta_xy = xy_size.getValue();
    inputs->LengthZ = thickness.getValue();
    inputs->tiltSelection = static_cast<SOC::TiltSelection>(tiltSelection.getValue());

    std::vector<uint8_t> viewMasks;

    if(viewMask.getValue().length() != 0 && parseUnknownArray(viewMask.getValue(), "%d", viewMasks) < 0)
    {
      std::cout << "Error Parsing the Tilt Mask Values. They should be entered as 4,7,12,34,67" << std::endl;
      return -1;
    }
    inputs->excludedViews = viewMasks;

  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << " error: " << e.error() << " for arg " << e.argId() << std::endl;
    std::cout << "** Unknown Arguments. Displaying help listing instead. **" << std::endl;
    return -1;
  }
  return 0;
}

#define PRINT_VAR(out, inputs, var)\
	out << #var << ": " << inputs->var << std::endl;
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCArgsParser::printInputs(TomoInputsPtr inputs, std::ostream &out)
{
  out << "------------------ TomoInputs Begin ------------------" << std::endl;
  PRINT_VAR(out, inputs, NumIter);
  PRINT_VAR(out, inputs, NumOuterIter);
  PRINT_VAR(out, inputs, SigmaX);
  PRINT_VAR(out, inputs, p);
  PRINT_VAR(out, inputs, StopThreshold);
  PRINT_VAR(out, inputs, InterpFlag);
  PRINT_VAR(out, inputs, useSubvolume);
  PRINT_VAR(out, inputs, xStart);
  PRINT_VAR(out, inputs, xEnd);
  PRINT_VAR(out, inputs, yStart);
  PRINT_VAR(out, inputs, yEnd);
  PRINT_VAR(out, inputs, zStart);
  PRINT_VAR(out, inputs, zEnd);
  PRINT_VAR(out, inputs, tiltSelection);
  PRINT_VAR(out, inputs, fileXSize);
  PRINT_VAR(out, inputs, fileYSize);
  PRINT_VAR(out, inputs, fileZSize);
  PRINT_VAR(out, inputs, LengthZ);
  PRINT_VAR(out, inputs, delta_xz);
  PRINT_VAR(out, inputs, delta_xy);
  PRINT_VAR(out, inputs, targetGain);
  PRINT_VAR(out, inputs, sinoFile);
  PRINT_VAR(out, inputs, initialReconFile);
  PRINT_VAR(out, inputs, gainsInputFile);
  PRINT_VAR(out, inputs, offsetsInputFile);
  PRINT_VAR(out, inputs, varianceInputFile);

  PRINT_VAR(out, inputs, tempDir);
  PRINT_VAR(out, inputs, reconstructedOutputFile);
  PRINT_VAR(out, inputs, gainsOutputFile);
  PRINT_VAR(out, inputs, offsetsOutputFile);
  PRINT_VAR(out, inputs, varianceOutputFile);

  out << "------------------ TomoInputs End ------------------" << std::endl;
}
