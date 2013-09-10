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



#include "MultiResolutionReconstruction.h"

#include <errno.h>

#include <iostream>

#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/MXAFileInfo.h"
#include "MXA/Utilities/StringUtils.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/HAADF/HAADFForwardModel.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionReconstruction::MultiResolutionReconstruction() :
FilterPipeline(),
m_Debug(false),
m_InputFile(""),
m_TempDir(""),
m_OutputFile(""),
m_BrightFieldFile(""),
m_InitialReconstructionFile(""),
m_DeleteTempFiles(false),
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
m_TiltSelection(SOC::A_Tilt),
m_Cancel(false)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MultiResolutionReconstruction::~MultiResolutionReconstruction()
{
}


#define PRINT_VAR(out, inputs, var)\
    out << #var << ": " << inputs->var << std::endl;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionReconstruction::printInputs(TomoInputsPtr inputs, std::ostream &out)
{
#if 1
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
  PRINT_VAR(out, inputs, extendObject);
  PRINT_VAR(out, inputs, interpolateFactor);
  PRINT_VAR(out, inputs, defaultInitialRecon);

#if 0
  PRINT_VAR(out, inputs, targetGain);
  PRINT_VAR(out, inputs, useDefaultOffset);
  PRINT_VAR(out, inputs, defaultOffset);
  PRINT_VAR(out, inputs, defaultVariance);
#endif

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
#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionReconstruction::setCancel(bool value)
{
 m_Cancel = value;
 if (NULL != m_CurrentEngine.get())
 {
   m_CurrentEngine->setCancel(value);
 }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool MultiResolutionReconstruction::getCancel()
{
  return m_Cancel;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MultiResolutionReconstruction::execute()
{
  std::vector<std::string>  tempFiles;

  pipelineProgressMessage("MultiResolutionSOC::execute");
  int err = 0;
  std::stringstream ss;
  ss << "-- There are " << m_NumberResolutions << " resolutions to reconstruct." ;
  pipelineProgressMessage(ss.str());
  ss.str("");

  TomoInputsPtr prevInputs = TomoInputsPtr(new TomoInputs);
  ReconstructionEngine::InitializeTomoInputs(prevInputs);


  for (int i = 0; i < m_NumberResolutions; ++i)
  {	  
    HAADFForwardModel::Pointer forwardModel = HAADFForwardModel::New();
    if(getCancel() == true)
    {
      setErrorCondition(-999);
      return;
    }

    TomoInputsPtr inputs = TomoInputsPtr(new TomoInputs);
    ReconstructionEngine::InitializeTomoInputs(inputs);

    /* ******* this is bad. Remove this for production work ****** */
    inputs->extendObject = getExtendObject();

    ss << "Extend Object Flag" << inputs->extendObject << std::endl;
    pipelineProgressMessage(ss.str());

    /* Get our input files from the last resolution iteration */
    inputs->gainsInputFile = prevInputs->gainsOutputFile;
    inputs->offsetsInputFile = prevInputs->offsetsOutputFile;
    inputs->varianceInputFile = prevInputs->varianceOutputFile;
    if(i == 0)
    {
      inputs->initialReconFile = getInitialReconstructionFile();
    }
    else
    {
      inputs->initialReconFile = prevInputs->reconstructedOutputFile;
    }

    if(i == 0)
    {
      inputs->InterpFlag = (getInterpolateInitialReconstruction() == false) ? 0 : 1;
    }
    else
    {
      inputs->InterpFlag = 1;
    } //If at a finer resolution need to interpolate

    inputs->interpolateFactor = powf((float)2, (float)getNumberResolutions() - 1) * m_FinalResolution;

    /* Now set the input files for this resolution */
    inputs->sinoFile = m_InputFile;
    inputs->tempDir = m_TempDir + MXADir::Separator + StringUtils::numToString(inputs->interpolateFactor / static_cast<int>(powf(2.0f, i))) + std::string("x");

    //Make sure the directory is created:
    bool success = MXADir::mkdir(inputs->tempDir, true);
    if(!success)
    {
      ss.str("");
      ss << "Could not create path: " << inputs->tempDir << std::endl;
      setErrorCondition(-1);
      pipelineErrorMessage(ss.str());
      return;
    }

	 //initialize the Bragg selector file
	  inputs->braggSelectorFile = MXAFileInfo::parentPath(m_OutputFile) + MXADir::Separator + MXAFileInfo::fileNameWithOutExtension(m_OutputFile) + "Selector.mrc";
	  
    // Only write the mrc and vtk files on the last iteration
    if(m_NumberResolutions - 1 == i) // Last Iteration
    {
        inputs->vtkOutputFile = MXAFileInfo::parentPath(m_OutputFile) + MXADir::Separator + MXAFileInfo::fileNameWithOutExtension(m_OutputFile) + ".vtk";
        inputs->avizoOutputFile = MXAFileInfo::parentPath(m_OutputFile) + MXADir::Separator + MXAFileInfo::fileNameWithOutExtension(m_OutputFile) + ".am";
        inputs->mrcOutputFile = m_OutputFile;
	
    }
    else
    {
        ss.str("");
        ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedObjectFile;
        inputs->reconstructedOutputFile = ss.str();

        inputs->vtkOutputFile = "";
        inputs->mrcOutputFile = "";
        inputs->avizoOutputFile = "";
    }

    // Line up all the temp files that we are going to delete
    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalGainParametersFile;
    inputs->gainsOutputFile = ss.str();
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalOffsetParametersFile;
    inputs->offsetsOutputFile = ss.str();
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::FinalVariancesFile;
    inputs->varianceOutputFile = ss.str();
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedSinogramFile;
    tempFiles.push_back(ss.str());

    // Create the paths for all the temp files that we want to delete
    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::CostFunctionFile;
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::DetectorResponseFile;
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::ReconstructedObjectFile;
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::VoxelProfileFile;
    tempFiles.push_back(ss.str());

    ss.str("");
    ss << inputs->tempDir << MXADir::Separator << ScaleOffsetCorrection::UpsampledBinFile;
    tempFiles.push_back(ss.str());


    inputs->NumOuterIter = getOuterIterations();
    if(i == 0)
    {

      inputs->NumIter = getInnerIterations();
    }
    else
    {
      inputs->NumIter = 1; //getInnerIterations();
    }
    inputs->p = getMRFShapeParameter() + 1; //NEW: Now we are going to compute p = Diffuseness(input) + 1
    inputs->StopThreshold = getStopThreshold();
    if(i >= 2)
    {
      inputs->StopThreshold = getStopThreshold() * 2.0f;
    }
    /** SIGMA_X needs to be calculated here based on some formula**/
    inputs->SigmaX = pow(2, (getNumberResolutions() - 1 - i) * (1 - 3 / inputs->p)) * getSigmaX();
    ss.str("");
    ss << "SigmaX=" << inputs->SigmaX;
    pipelineProgressMessage(ss.str());
    if(i == 0)
    {
      inputs->defaultInitialRecon = getInitialReconstructionValue();
      forwardModel->setDefaultVariance(getDefaultVariance());
    }
    inputs->delta_xy = powf(2.0f, getNumberResolutions() - i - 1) * static_cast<Real_t>(m_FinalResolution);
    inputs->delta_xz = powf(2.0f, getNumberResolutions() - i - 1) * static_cast<Real_t>(m_FinalResolution);

    if(i == 0)
    {
      forwardModel->setDefaultOffset(getDefaultOffsetValue());
      forwardModel->setUseDefaultOffset(getUseDefaultOffset());
    }
    inputs->LengthZ = m_SampleThickness;
    forwardModel->setTargetGain(getTargetGain());//In the BF code this is always set to 1
	
	forwardModel->setBraggDelta(getBraggDelta()); //Set the Bragg function Delta value 
    
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

    }
	  

    inputs->excludedViews = m_ViewMasks;

    SinogramPtr sinogram = SinogramPtr(new Sinogram);
    GeometryPtr geometry = GeometryPtr(new Geometry);
  

    ReconstructionEngine::InitializeSinogram(sinogram);
    ReconstructionEngine::InitializeGeometry(geometry);

    //Calculate approximate memory required
	  memCalculate(inputs);

    forwardModel->setAdvParams(m_AdvParams);
    forwardModel->setTomoInputs(inputs);

    forwardModel->addObserver(this);
    forwardModel->setMessagePrefix(StringUtils::numToString(inputs->interpolateFactor / static_cast<int>(powf(2.0f, i))) + std::string("x: "));

    forwardModel->setVerbose(true);
    forwardModel->setVeryVerbose(true);

    //Create an Engine and initialize all the structures
    ReconstructionEngine::Pointer engine = ReconstructionEngine::New();
    m_CurrentEngine = engine;
    engine->setTomoInputs(inputs);
    engine->setSinogram(sinogram);
    engine->setGeometry(geometry);
    engine->setAdvParams(m_AdvParams);
    engine->setForwardModel(forwardModel);

    engine->setVerbose(true);
    engine->setVeryVerbose(true);


    // We need to get messages to the gui or command line
    engine->addObserver(this);
    engine->setMessagePrefix(StringUtils::numToString(inputs->interpolateFactor / static_cast<int>(powf(2.0f, i))) + std::string("x: "));
    ss.str("");
    ss << "Sinogram Inputs -----------------------------------------" << std::endl;
    printInputs(inputs, ss);
    pipelineProgressMessage(ss.str());
    
    engine->execute();
    engine = ReconstructionEngine::NullPointer();

    prevInputs = inputs;

    // Get any tempfiles created by the process such as intermediate files for display
    // during the reconstruction
    for(size_t i = 0; i < inputs->tempFiles.size(); ++i)
    {
        tempFiles.push_back(inputs->tempFiles[i]);
    }

    // Tack on the output temp directory last since we delete in the order we placed them
    // in the vector. This will make sure the directory is empty before we try to delete it.
    // If the directory is NOT empty then the user added something to it that we don't know
    // about so DO NOT try to delete the directory
    ss.str("");
    ss << inputs->tempDir;
    tempFiles.push_back(ss.str());
  }


  if (getDeleteTempFiles() == true)
  {
    for(size_t i = 0; i < tempFiles.size(); ++i)
    {
       // std::cout << "Removing: " << tempFiles[i] << std::endl;
        if(MXADir::isDirectory(tempFiles[i]) == true )
        {
            errno = 0;
            if (false == MXADir::rmdir(tempFiles[i], false) )
            {
                std::cout << errno << " - Could NOT remove Directory: " << tempFiles[i] << std::endl;
                std::vector<std::string> dirList = MXADir::entryList(tempFiles[i]);
                for(size_t i = 0; i < dirList.size(); ++i)
                {
                    std::cout << "   " << dirList[i] << std::endl;
                }
            }
        }
        else
        {
            MXAFileInfo::remove(tempFiles[i]);
        }
    }
  }

  updateProgressAndMessage("MultiResolution SOC Complete", 100);
  setErrorCondition(err);
}

void MultiResolutionReconstruction::memCalculate(TomoInputsPtr inputs)
{
    float GeomNx,GeomNy,GeomNz;
    float SinoNr,SinoNt,SinoNtheta;
    SinoNr = static_cast<float>(inputs->xEnd - inputs->xStart+1);
    SinoNt = static_cast<float>(inputs->yEnd - inputs->yStart+1);
    SinoNtheta = static_cast<float>(inputs->zEnd - inputs->zStart+1);

    AdvancedParametersPtr advancedParams = AdvancedParametersPtr(new AdvancedParameters);
    ReconstructionEngine::InitializeAdvancedParams(advancedParams);

    if(inputs->extendObject == 1)
    {
        GeomNx = (SinoNr/m_FinalResolution)*4;//TODO:Need to access X_Stretch and
    }
    else
    {
        GeomNx = SinoNr/m_FinalResolution;
    }

    GeomNy = SinoNt/m_FinalResolution;
    GeomNz = advancedParams->Z_STRETCH*(m_SampleThickness/(m_FinalResolution));// TODO: need to access Sinogram_deltar and z_stretch.
    //This is wrong currently. Need to multiply m_FinalResolution by size of voxel in nm

    float dataTypeMem = sizeof(Real_t);
    float ObjectMem = GeomNx*GeomNy*GeomNz*dataTypeMem;
    float SinogramMem = SinoNr*SinoNt*SinoNtheta*dataTypeMem;
    float ErroSinoMem = SinogramMem;
    float WeightMem = SinogramMem; //Weight matrix
    float A_MatrixMem;
    if(0 == inputs->extendObject)
    {
    A_MatrixMem = GeomNx*GeomNz*(m_FinalResolution*3*(dataTypeMem+4)*SinoNtheta);// 4 is the bytes to store the counts
   //*+4 correspodns to bytes to store a single double and a unsigned into to
   //store the offset. 3*m_FinalRes is the approximate number of detector elements hit per voxel
    }
    else {
        A_MatrixMem = GeomNx*GeomNz*(m_FinalResolution*(dataTypeMem+4)*SinoNtheta); //Since we are reconstructing a larger region there are several voxels with no projection data. so instead of each voxel hitting 3*m_FinalRes det entries we aproximate it by m_FinalRes
    }
    float NuisanceParamMem = SinoNtheta*dataTypeMem*3;//3 is for gains offsets and noise var

    float TotalMem = ObjectMem+SinogramMem*2+ErroSinoMem+WeightMem+A_MatrixMem+NuisanceParamMem;//in bytes

    TotalMem/=(1e9);//To get answer in Gb

    std::cout<<"Total Max Mem needed = "<<TotalMem<<" Gb"<<std::endl;

}
