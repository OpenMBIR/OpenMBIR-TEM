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
#include "BFReconstructionArgsParser.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>

#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/MXAFileInfo.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/MBIRLibVersion.h"
#include "MBIRLib/Common/EIMTime.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Common/allocate.h"
#include "MBIRLib/IOFilters/MRCHeader.h"
#include "MBIRLib/IOFilters/MRCReader.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/BrightField/BFConstants.h"


/**
 * @brief Parses numeric values from a delimited string into a preallocated array storage.
 * The programmer MUST know in advance how many values there will be.
 * @param values The string to be parsed
 * @param format The stdio format specifier to use (%f for floats, %d for integers
 * @param output The output location to store the parsed values
 * @return Error condition
 */
template<typename T>
int parseValues(const std::string& values, const char* format, T* output)
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
int parseUnknownArray(const std::string& values, const char* format, std::vector<T>& output)
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
BFReconstructionArgsParser::BFReconstructionArgsParser()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BFReconstructionArgsParser::~BFReconstructionArgsParser()
{

}


// -----------------------------------------------------------------------------
// Parse the input arguments in the command line tool
// -----------------------------------------------------------------------------
int BFReconstructionArgsParser::parseArguments(int argc, char** argv, BFMultiResolutionReconstruction::Pointer m_MultiResSOC)
{
  if(NULL == m_MultiResSOC)
  {
    printf("The MultiResolutionSOC was null. Returning early.\n");
    return -1;
  }

  TCLAP::CmdLine cmd("", ' ', MBIRLib::Version::Complete());

  TCLAP::ValueArg<std::string> in_MRCFile("s", "sinofile", "The Sinogram File", true, "", "");
  cmd.add(in_MRCFile);

  // TCLAP::ValueArg<std::string> inputBrightFieldFilePath("", "brightfield", "Acquired Bright Field tilt series", false, "", "");
  // cmd.add(inputBrightFieldFilePath);

  TCLAP::ValueArg<std::string> reconstructedVolumeFileName("", "outputfile", "The Output File", true, "", "");
  cmd.add(reconstructedVolumeFileName);


  TCLAP::ValueArg<std::string> subvolume("", "subvolume", "SubVolume to Reconstruct in the form xmin,ymin,zmin,xmax,ymax,zmax", true, "", "");
  cmd.add(subvolume);
  TCLAP::ValueArg<Real_t> sampleThickness("", "thickness", "Thickness of sample in nm", true, 0, "");
  cmd.add(sampleThickness);
  TCLAP::ValueArg<unsigned int> tiltSelection("", "tilt_selection", "Which Tilt Values to use from file. Default is 'A'", false, SOC::A_Tilt, "0");
  cmd.add(tiltSelection);
  TCLAP::ValueArg<double> mrf("", "diffuseness", "Diffuseness Parameter (Markov Random Field Parameter)", false, 0.200, "0.2");
  cmd.add(mrf);
  TCLAP::ValueArg<int> finalResolution("", "final_resolution", "The final resolution multiple for the reconstruction" , false, 1, "1");
  cmd.add(finalResolution);
  TCLAP::ValueArg<double> sigma_x("", "sigma_x", "Sigma X Value", true, 1.0, "1.0");
  cmd.add(sigma_x);


  TCLAP::ValueArg<double> stopThreshold("", "stop_threshold", "Stopping Threshold for inner loop", false, .001, ".001");
  cmd.add(stopThreshold);
  TCLAP::ValueArg<int> numResolutions("", "num_resolutions", "The number of resolutions to use" , false, 3, "3");
  cmd.add(numResolutions);
  TCLAP::ValueArg<double> targetGain("", "target_gain", "Target Gain for unscattered electrons", false, 1, "1");
  cmd.add(targetGain);

  TCLAP::ValueArg<double> braggThreshold("", "bragg_T", "Value of Bragg threshold T", false, 3, "3");
  cmd.add(braggThreshold);

  TCLAP::ValueArg<double> braggDelta("", "bragg_delta", "value of delta in the bragg rejector", false, 0.5, "0.5");
  cmd.add(braggDelta);

  TCLAP::ValueArg<double> bfOffset("", "bf_offset", "value of offset in the BF data", false, 0.0, "0");
  cmd.add(bfOffset);

  TCLAP::SwitchArg interpolateInitialRecontruction ("", "interpolate_initial_recon", "Interpolate Initial Reconstruction Value", false);
  cmd.add(interpolateInitialRecontruction);
  TCLAP::ValueArg<int> outerIterations("", "outer_iterations", "Outer Iterations to use", false, 600, "600");
  cmd.add(outerIterations);
  TCLAP::ValueArg<int> innerIterations("", "inner_iterations", "Number of Inner Iterations", false, 100, "100");
  cmd.add(innerIterations);
  TCLAP::ValueArg<double> defaultOffset("", "default_dosage", "Default dosage for all Tilts", true, 0, "");
  cmd.add(defaultOffset);
  TCLAP::SwitchArg useDefaultOffset("", "use_default_dosage", "Use the Default dosage Value" , false);
  cmd.add(useDefaultOffset);
  TCLAP::ValueArg<double> defaultVariance("", "default_variance", "Default variance for all Tilts", false, 1.0, "1");
  cmd.add(defaultVariance);
  TCLAP::ValueArg<double> defaultInitialRecon("", "default_recon_value", "Default initial value of reconstruction", false, 0.0, "0");
  cmd.add(defaultInitialRecon);
  TCLAP::SwitchArg extendObject("", "extend_object", "To extend the object or not", false);
  cmd.add(extendObject);
  TCLAP::SwitchArg m_DeleteTempFiles ("", "delete_tmp_files", "Delete all the Temp files that are created", false);
  cmd.add(m_DeleteTempFiles);


  TCLAP::ValueArg<std::string> initialReconstructionPath("i", "initial_recon_file", "Initial Reconstruction to initialize algorithm", false, "", "");
  cmd.add(initialReconstructionPath);

#if 0
  TCLAP::ValueArg<std::string> in_Gains("", "gains", "Initial Gains to use.", false, "", "");
  cmd.add(in_Gains);

  TCLAP::ValueArg<std::string> in_Offsets("", "offsets", "Initial Offsets to use.", false, "", "");
  cmd.add(in_Offsets);

  TCLAP::ValueArg<std::string> in_Variance("", "variance", "Initial Variance to use.", false, "", "");
  cmd.add(in_Variance);

  TCLAP::ValueArg<double> in_InterpFactor("", "interpFactor", "Interpolate Factor", false, 0.0, "");
  cmd.add(in_InterpFactor);

  TCLAP::ValueArg<Real_t> xz_size("", "xz_size", "Size in nm of output pixel xz plane", true, 1, "1");
  cmd.add(xz_size);
  TCLAP::ValueArg<Real_t> xy_size("", "xy_size", "Size in nm of output pixel xy plane", true, 1, "1");
  cmd.add(xy_size);
#endif
  TCLAP::ValueArg<std::string> viewMask("", "exclude_views", "Comma separated list of tilts to exclude by index", false, "", "");
  cmd.add(viewMask);


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

    std::string path;
    path = MXADir::toNativeSeparators(in_MRCFile.getValue());
    m_MultiResSOC->setInputFile(path);

    path = MXAFileInfo::parentPath(reconstructedVolumeFileName.getValue());
    m_MultiResSOC->setTempDir(path);

    path = MXADir::toNativeSeparators(MXAFileInfo::absolutePath(reconstructedVolumeFileName.getValue()));
    m_MultiResSOC->setOutputFile(path);

    path = MXADir::toNativeSeparators(initialReconstructionPath.getValue());
    m_MultiResSOC->setInitialReconstructionFile(path);

    m_MultiResSOC->setNumberResolutions(numResolutions.getValue());

    m_MultiResSOC->setFinalResolution(finalResolution.getValue());
    m_MultiResSOC->setSampleThickness(sampleThickness.getValue());
    m_MultiResSOC->setTargetGain(targetGain.getValue());
    m_MultiResSOC->setBraggThreshold(braggThreshold.getValue());
    m_MultiResSOC->setBraggDelta(braggDelta.getValue());
    m_MultiResSOC->setBfOffset(bfOffset.getValue());
    m_MultiResSOC->setStopThreshold(stopThreshold.getValue());
    m_MultiResSOC->setOuterIterations(outerIterations.getValue());
    m_MultiResSOC->setInnerIterations(innerIterations.getValue());
    m_MultiResSOC->setSigmaX(sigma_x.getValue());

    m_MultiResSOC->setMRFShapeParameter(mrf.getValue());

    //Based on the dosage in a void and the offset in the data compute the "d" parameter
    Real_t true_offset = -log(defaultOffset.getValue() + bfOffset.getValue());
    m_MultiResSOC->setDefaultOffsetValue(true_offset);

    m_MultiResSOC->setUseDefaultOffset(useDefaultOffset.getValue());
    //    m_MultiResSOC->setExtendObject(extendObject.getValue());
    m_MultiResSOC->setDefaultVariance(defaultVariance.getValue());
    m_MultiResSOC->setInitialReconstructionValue(defaultInitialRecon.getValue());

    m_MultiResSOC->setInterpolateInitialReconstruction(interpolateInitialRecontruction.getValue());
    m_MultiResSOC->setDeleteTempFiles(m_DeleteTempFiles.getValue());
    AdvancedParametersPtr advParams = AdvancedParametersPtr(new AdvancedParameters);
    BFReconstructionEngine::InitializeAdvancedParams(advParams);
    m_MultiResSOC->setAdvParams(advParams);

    int subvolumeValues[6];
    ::memset(subvolumeValues, 0, 6 * sizeof(int));

    if(subvolume.getValue().length() != 0)
    {
      int err = parseValues(subvolume.getValue(), "%d", subvolumeValues);
      if(err < 0)
      {
        std::cout << "Error Parsing the Subvolume Dimensions. They should be entered as --subvolume 64,128,256,80,150,280" << std::endl;
        return -1;
      }
    }

    std::vector<uint16_t> multiResSubVolume(6);
    for(int i = 0; i < 6; ++i) { multiResSubVolume[i] = subvolumeValues[i]; }
    m_MultiResSOC->setSubvolume(multiResSubVolume);


    std::vector<uint8_t> viewMasks;

    if(viewMask.getValue().length() != 0 && parseUnknownArray(viewMask.getValue(), "%d", viewMasks) < 0)
    {
      std::cout << "Error Parsing the Tilt Mask Values. They should be entered as 4,7,12,34,67" << std::endl;
      return -1;
    }

    m_MultiResSOC->setViewMasks(viewMasks);



    // Read the proper tilts from the mrc file
    MRCHeader header;
    ::memset(&header, 0, sizeof(header));
    MRCReader::Pointer reader = MRCReader::New(true);
    // Read the header from the file
    int err = reader->readHeader(m_MultiResSOC->getInputFile(), &header);
    if(err < 0)
    {
      return -1;
    }
    //  int tiltIndex = 0;
    if(header.feiHeaders != NULL)
    {
      std::vector<float> tilts(header.nz, 0.0f);
      //FEIHeader fei = header.feiHeaders[tiltIndex];
      for(int l = 0; l < header.nz; ++l)
      {
        if(tiltSelection.getValue() == 0)
        {
          tilts[l] = header.feiHeaders[l].a_tilt;
        }
        else
        {
          tilts[l] = header.feiHeaders[l].b_tilt;
        }
      }
      m_MultiResSOC->setTilts(tilts);
    }

  }
  catch (TCLAP::ArgException& e)
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
void BFReconstructionArgsParser::printInputs(TomoInputsPtr inputs, std::ostream& out)
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
  //PRINT_VAR(out, inputs, tiltSelection);
  PRINT_VAR(out, inputs, fileXSize);
  PRINT_VAR(out, inputs, fileYSize);
  PRINT_VAR(out, inputs, fileZSize);
  PRINT_VAR(out, inputs, LengthZ);
  PRINT_VAR(out, inputs, delta_xz);
  PRINT_VAR(out, inputs, delta_xy);
  PRINT_VAR(out, inputs, extendObject);
  PRINT_VAR(out, inputs, interpolateFactor);
  PRINT_VAR(out, inputs, defaultInitialRecon);
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
