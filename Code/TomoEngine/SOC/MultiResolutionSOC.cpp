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



#include "MultiResolutionSOC.h"

#include <iostream>

#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/StringUtils.h"
#include "TomoEngine/Common/EIMMath.h"

#include "TomoEngine/SOC/SOCEngine.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionSOC::MultiResolutionSOC() :
m_InputFile(""),
m_TempDir(""),
m_OutputFile(""),
m_BrightFieldFile(""),
m_InitialReconstructionFile(""),
m_NumberResolutions(1),
m_SampleThickness(100.0f),
m_TargetGain(0.0f),
m_StopThreshold(0.009f),
m_OuterIterations(1),
m_InnerIterations(1),
m_SigmaX(0.0f),
m_MRFShapeParameter(1.1f),
m_DefaultOffsetValue(0.0f),
m_UseDefaultOffset(false),
m_FinalResolution(1),
m_ExtendObject(true),
m_InterpolateInitialReconstruction(false),
m_DefaultVariance(1.0f),
m_InitialReconstructionValue(0.0f),
m_TiltSelection(SOC::A_Tilt)

{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionSOC::~MultiResolutionSOC()
{
}


#define PRINT_VAR(out, inputs, var)\
	out << #var << ": " << inputs->var << std::endl;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionSOC::printInputs(TomoInputsPtr inputs, std::ostream &out)
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

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionSOC::execute()
{
  std::cout << "MultiResolutionSOC::execute" << std::endl;
  int err = 0;
  std::cout << "-- There are " << m_NumberResolutions << " resolutions to reconstruct." << std::endl;

  std::stringstream ss;


  TomoInputsPtr prevInputs = TomoInputsPtr(new TomoInputs);
  SOCEngine::InitializeTomoInputs(prevInputs);

  TomoInputsPtr bf_inputs = TomoInputsPtr(new TomoInputs);
  SOCEngine::InitializeTomoInputs(bf_inputs);

  for (int i = 0; i < m_NumberResolutions; ++i)
  {
    TomoInputsPtr inputs = TomoInputsPtr(new TomoInputs);
    SOCEngine::InitializeTomoInputs(inputs);

    bf_inputs->sinoFile = getBrightFieldFile();

    /* ******* this is bad. Remove this for production work ****** */
    inputs->extendObject = getExtendObject();

	  std::cout<<"Extend Object Flag"<<inputs->extendObject<<std::endl;

    /* Get our input files from the last resolution iteration */
    inputs->gainsInputFile = prevInputs->gainsOutputFile;
    inputs->offsetsInputFile = prevInputs->offsetsOutputFile;
    inputs->varianceInputFile = prevInputs->varianceOutputFile;
    if ( i == 0) {
      inputs->initialReconFile = getInitialReconstructionFile();
    }
    else
    {
      inputs->initialReconFile = prevInputs->reconstructedOutputFile;
    }

    if (i == 0)
    {
      inputs->InterpFlag = (getInterpolateInitialReconstruction() == false) ? 0 : 1;
    }
    else
    { inputs->InterpFlag = 1; }//If at a finer resolution need to interpolate


    inputs->interpolateFactor = powf((float)2, (float)getNumberResolutions()-1) * m_FinalResolution;

    /* Now set the output files for this resolution */
    inputs->sinoFile = m_InputFile;
    inputs->tempDir  = m_TempDir + MXADir::Separator + StringUtils::numToString(inputs->interpolateFactor/(powf(2.0f,i))) + std::string("x");

    //Make sure the directory is created:
    bool success = MXADir::mkdir(inputs->tempDir, true);
    if (!success)
    {
      std::cout << "Could not create path: " << inputs->tempDir << std::endl;
    }

	  /***************TO DO - Fix This*********************/
    if(m_NumberResolutions - 1 == i)
    {
      ss.str("");
      ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedBinFile;
      inputs->reconstructedOutputFile = ss.str();
    }
    else
    {
      ss.str("");
      ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedBinFile;
      inputs->reconstructedOutputFile = ss.str();
    }
    /************************************/

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalGainParametersFile;
    inputs->gainsOutputFile = ss.str();

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalOffsetParametersFile;
    inputs->offsetsOutputFile = ss.str();

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalVariancesFile;
    inputs->varianceOutputFile = ss.str();

    inputs->NumOuterIter = getOuterIterations();
    if(i == 0)
    {
//      inputs->NumIter = 20;
      inputs->NumIter = getInnerIterations();
    }
    else
    {
		inputs->NumIter = 1;//getInnerIterations();
    }
    inputs->p = getMRFShapeParameter();
    inputs->StopThreshold = getStopThreshold();
    if (i >= 2)
    {
      inputs->StopThreshold = getStopThreshold()*2.0f;
    }
    /** SIGMA_X needs to be calculated here based on some formula**/
    inputs->SigmaX =pow(2,(getNumberResolutions()-1-i)*(1-3/inputs->p)) * getSigmaX();
	  std::cout<<"SigmaX="<<inputs->SigmaX;
	  if (i == 0)
	  {
	    inputs->defaultInitialRecon = getInitialReconstructionValue();
      inputs->defaultVariance = getDefaultVariance();
	  }
    inputs->delta_xy = powf(2.0f, getNumberResolutions()-i-1)*m_FinalResolution;
    inputs->delta_xz = powf(2.0f, getNumberResolutions()-i-1)*m_FinalResolution;

    if (i == 0)
    {
      inputs->defaultOffset = getDefaultOffsetValue();
      inputs->useDefaultOffset = getUseDefaultOffset();
    }
    inputs->LengthZ = m_SampleThickness;
    inputs->targetGain = m_TargetGain;
    inputs->tiltSelection = m_TiltSelection;
    if(m_Subvolume.size() > 0)
    {
      inputs->useSubvolume = true;
      inputs->xStart = m_Subvolume[0];
      inputs->xEnd = m_Subvolume[3];
      inputs->yStart = m_Subvolume[1];
      inputs->yEnd = m_Subvolume[4];
      inputs->zStart = m_Subvolume[2];
      inputs->zEnd = m_Subvolume[5];

      bf_inputs->useSubvolume = true;
      bf_inputs->xStart = m_Subvolume[0];
      bf_inputs->xEnd = m_Subvolume[3];
      bf_inputs->yStart = m_Subvolume[1];
      bf_inputs->yEnd = m_Subvolume[4];
      bf_inputs->zStart = m_Subvolume[2];
      bf_inputs->zEnd = m_Subvolume[5];
    }
    inputs->excludedViews = m_ViewMasks;
    bf_inputs->excludedViews = m_ViewMasks;

    //Adjusting the volume along the y-directions so we dont have
    //  issues with pixelation
    int16_t disty = inputs->yEnd - inputs->yStart + 1;
    //3*iterpFactor is to account for the prior which operates on
    //26 point 3-D neighborhood which needs 3 x-z slices at the least
    int16_t rem_temp = disty % ((int16_t)inputs->interpolateFactor * 3);
    if(rem_temp != 0)
    {
      std::cout << "The number of y-pixels is not a proper multiple for multi-res" << std::endl;
      int16_t remainder = static_cast<int16_t>((inputs->interpolateFactor * 3) - (rem_temp));
      inputs->yEnd += remainder;
    }

    SinogramPtr sinogram = SinogramPtr(new Sinogram);
    SinogramPtr bf_sinogram = SinogramPtr(new Sinogram);
    GeometryPtr geometry = GeometryPtr(new Geometry);
    ScaleOffsetParamsPtr nuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);

    SOCEngine::InitializeSinogram(sinogram);
    SOCEngine::InitializeGeometry(geometry);
    SOCEngine::InitializeScaleOffsetParams(nuisanceParams);
    SOCEngine::InitializeSinogram(bf_sinogram);

    // Create an Engine and initialize all the structures
	  SOCEngine::Pointer engine = SOCEngine::New();
	  engine->setTomoInputs(inputs);
    engine->setSinogram(sinogram);
    engine->setGeometry(geometry);
    engine->setNuisanceParams(nuisanceParams);
    engine->setBFTomoInputs(bf_inputs);
    engine->setBFSinogram(bf_sinogram);
    // We need to get messages to the gui or command line
    engine->addObserver(this);
    engine->setMessagePrefix( StringUtils::numToString(inputs->interpolateFactor/(powf(2.0f,i))) + std::string("x: ") );

    printInputs(inputs, std::cout);
	  
	printInputs(bf_inputs, std::cout);  

    engine->execute();
    engine = SOCEngine::NullPointer();

    prevInputs = inputs;
  }

  updateProgressAndMessage("MultiResolution SOC Complete", 100);
  setErrorCondition(err);
}
