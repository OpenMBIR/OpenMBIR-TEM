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
#include "TomoEngine/Common/EIMTime.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/IO/MRCReader.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/TomoEngineVersion.h"



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
  int n = sscanf(values.substr(0, pos).c_str(), format, &(output[index]) );
  if (n != 1)
  {
    return -1;
  }

  ++index;
  while(pos != std::string::npos && pos != values.size() - 1)
  {
    n = sscanf(values.substr(pos+1).c_str(), format, &(output[index]) );
    pos = values.find(",", pos+1);
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
  int n = sscanf(values.substr(0, pos).c_str(), format, &t );
  if (n != 1)
  {
    return -1;
  }
  output.push_back(t);

  ++index;
  while(pos != std::string::npos && pos != values.size() - 1)
  {
    n = sscanf(values.substr(pos+1).c_str(), format, &(t) );
    output.push_back(t);
    pos = values.find(",", pos+1);
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



int SOCArgsParser::parseArguments(int argc,char **argv, TomoInputs* Input)
{
  if ( NULL == Input)
  {
    printf("The CommandLineInputspointer was null. Returning early.\n");
    return -1;
  }

  TCLAP::CmdLine cmd("", ' ', TomoEngine::Version::Complete);

  TCLAP::ValueArg<std::string> in_sinoFile("s", "sinofile", "The Sinogram File", true, "", "");
  cmd.add(in_sinoFile);

  TCLAP::ValueArg<std::string> in_InitialRecon("i", "ini_recon", "Initial Reconstruction to initialize algorithm", false, "", "");
  cmd.add(in_InitialRecon);

  TCLAP::ValueArg<std::string> in_GainsOffsets("", "gains_offsets", "Initial Gains and Offsets to use.", false, "", "");
  cmd.add(in_GainsOffsets);

  TCLAP::ValueArg<std::string> in_outputFile("o", "outputfile", "The Output File", true, "", "");
  cmd.add(in_outputFile);

  TCLAP::ValueArg<std::string> in_outputDir("", "outdir", "The Output dir", true, ".", ".");
  cmd.add(in_outputDir);

	TCLAP::ValueArg<double> in_stopThreshold("T","stopThreshold","Stopping Threshold for inner loop",false, .009, "");
	cmd.add(in_stopThreshold);


  TCLAP::ValueArg<std::string> subvolume("", "subvolume", "SubVolume to Reconstruct in the form xmin,ymin,zmin,xmax,ymax,zmax", false, "", "");
  cmd.add(subvolume);



  TCLAP::ValueArg<int> in_numIter("n", "numIter", "Number of Iterations", true, 0, "0");
  cmd.add(in_numIter);
  TCLAP::ValueArg<double> in_sigmaX("l", "sigmax", "Sigma X Value", true, 1.0, "1.0");
  cmd.add(in_sigmaX);

  TCLAP::ValueArg<double> in_markov("m", "mrf", "Markov Random Field Parameter", true, 0.0, "0.0");
  cmd.add(in_markov);

	//TCLAP::ValueArg<std::string> InitialParameters("g", "initp", "InitialParameters", false, "", "");
 // cmd.add(InitialParameters);

  TCLAP::ValueArg<int> NumOuterIter("O", "", "NumOuterIter", true, 1, "1");
  cmd.add(NumOuterIter);

  TCLAP::ValueArg<DATA_TYPE> xz_size("", "xz_size", "Size in nm of output pixel xz plane", true, 1, "1");
  cmd.add(xz_size);
  TCLAP::ValueArg<DATA_TYPE> xy_size("", "xy_size", "Size in nm of output pixel xy plane", true, 1, "1");
  cmd.add(xy_size);

  TCLAP::ValueArg<DATA_TYPE> thickness("", "thickness", "Thickness of sample in nm", true, 0, "");
  cmd.add(thickness);

  TCLAP::ValueArg<std::string> viewMask("", "exclude_views", "Comma separated list of tilts to exclude by index", false, "", "");
  cmd.add(viewMask);

  if (argc < 2)
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
    Input->OutputFile = in_outputDir.getValue() + MXADir::getSeparator() + in_outputFile.getValue();
    Input->SinoFile = in_sinoFile.getValue();
    Input->InitialReconFile = in_InitialRecon.getValue();
    Input->GainsOffsetsFile = in_GainsOffsets.getValue();
    Input->outputDir = in_outputDir.getValue();
    Input->NumIter = in_numIter.getValue();
    Input->SigmaX = in_sigmaX.getValue();
    Input->p = in_markov.getValue();
    Input->NumOuterIter = NumOuterIter.getValue();
    Input->StopThreshold = in_stopThreshold.getValue();

    int subvolumeValues[6];
    ::memset(subvolumeValues, 0, 6 * sizeof(int));
    Input->useSubvolume = false;

    if(subvolume.getValue().length() != 0)
    {
      int err = parseValues(subvolume.getValue(), "%d", subvolumeValues);
      if(err < 0)
      {
        std::cout << "Error Parsing the Subvolume Dimensions. They should be entered as --subvolume 64,128,256,80,150,280" << std::endl;
        return -1;
      }
      Input->useSubvolume = true;
      Input->xStart = subvolumeValues[0];
      Input->xEnd = subvolumeValues[3];
      Input->yStart = subvolumeValues[1];
      Input->yEnd = subvolumeValues[4];
      Input->zStart = subvolumeValues[2];
      Input->zEnd = subvolumeValues[5];
    }

    Input->delta_xz = xz_size.getValue();
    Input->delta_xy = xy_size.getValue();
    Input->LengthZ = thickness.getValue();

    std::vector<uint8_t> viewMasks;

    if( viewMask.getValue().length() != 0
         && parseUnknownArray(viewMask.getValue(), "%d", viewMasks) < 0)
    {
      std::cout << "Error Parsing the Tilt Mask Values. They should be entered as 4,7,12,34,67" << std::endl;
      return -1;
    }
    Input->excludedViews = viewMasks;

  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << " error: " << e.error() << " for arg " << e.argId() << std::endl;
    std::cout << "** Unknown Arguments. Displaying help listing instead. **" << std::endl;
    return -1;
  }
  return 0;
}


#define PRINT_INPUT_VAR(input, var)\
  std::cout << #var << ": " << input->var << std::endl;
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCArgsParser::printArgs(std::ostream &out, TomoInputs* Input)
{
  out << "SOC Inputs ------------------------" << std::endl;
  PRINT_INPUT_VAR(Input, SinoFile)
  PRINT_INPUT_VAR(Input, InitialReconFile)
  PRINT_INPUT_VAR(Input, GainsOffsetsFile)
  PRINT_INPUT_VAR(Input, OutputFile)
  PRINT_INPUT_VAR(Input, outputDir)
  PRINT_INPUT_VAR(Input, NumIter)
  PRINT_INPUT_VAR(Input, NumOuterIter)
  PRINT_INPUT_VAR(Input, SigmaX)
  PRINT_INPUT_VAR(Input, p)
  PRINT_INPUT_VAR(Input, StopThreshold)
  PRINT_INPUT_VAR(Input, useSubvolume)
  PRINT_INPUT_VAR(Input, xStart)
  PRINT_INPUT_VAR(Input, xEnd)
  PRINT_INPUT_VAR(Input, yStart)
  PRINT_INPUT_VAR(Input, yEnd)
  PRINT_INPUT_VAR(Input, zStart)
  PRINT_INPUT_VAR(Input, zEnd)
  PRINT_INPUT_VAR(Input, LengthZ)
  PRINT_INPUT_VAR(Input, delta_xz)
  PRINT_INPUT_VAR(Input, delta_xy)

}
