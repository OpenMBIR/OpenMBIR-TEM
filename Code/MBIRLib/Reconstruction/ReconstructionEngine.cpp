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
#include "ReconstructionEngine.h"

// C Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// C++ Includes
#include <limits>
#include <iostream>


// MXA includes
#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/MXAFileInfo.h"

// Our own includes
#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Common/allocate.h"
#include "MBIRLib/Common/EIMTime.h"
#include "MBIRLib/Common/CE_ConstraintEquation.hpp"
#include "MBIRLib/Common/DerivOfCostFunc.hpp"

#include "MBIRLib/Reconstruction/ReconstructionConstants.h"

#include "MBIRLib/IOFilters/MRCHeader.h"
#include "MBIRLib/IOFilters/MRCReader.h"
#include "MBIRLib/IOFilters/MRCWriter.h"
#include "MBIRLib/IOFilters/RawGeometryWriter.h"

#include "MBIRLib/IOFilters/VTKFileWriters.hpp"
#include "MBIRLib/IOFilters/AvizoUniformCoordinateWriter.h"
#include "MBIRLib/IOFilters/DetectorResponseWriter.h"

#include "MBIRLib/GenericFilters/TomoFilter.h"
#include "MBIRLib/GenericFilters/DetectorResponse.h"

#include "MBIRLib/GenericFilters/RawSinogramInitializer.h"
#include "MBIRLib/GenericFilters/InitialReconstructionInitializer.h"
#include "MBIRLib/GenericFilters/InitialReconstructionBinReader.h"




#define USE_TBB_TASK_GROUP 1
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#include <tbb/task.h>
#endif

#define START_TIMER startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME(msg, startm, stopm)\
    std::cout << indent << msg << ": " << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

#define MAKE_OUTPUT_FILE(Fp, outdir, filename)\
    {\
    std::string filepath(outdir);\
    filepath = filepath.append(MXADir::getSeparator()).append(filename);\
    errno = 0;\
    Fp = fopen(filepath.c_str(),"wb");\
    if (Fp == NULL || errno > 0) { std::cout << "Error " << errno << " Opening Output file " << filepath << std::endl;}\
    }

#define COPY_333_ARRAY(i_max, j_max, k_max, src, dest)\
for(int i = 0; i < i_max; ++i){\
for(int j = 0; j < j_max; ++j){\
for(int k = 0; k < k_max; ++k){\
  dest[i][j][k] = src[i][j][k];\
}}}

namespace Detail {
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  inline double Clip(double x, double a, double b)
  {
    return (x < a) ? a : ((x > b) ? b:x);
  }


  // -----------------------------------------------------------------------------
  //<
  // -----------------------------------------------------------------------------
  inline int16_t mod(int16_t a,int16_t b)
  {
    int16_t temp;
    temp=a%b;
    if(temp < 0)
      return temp + b;
    else {
      return temp;
    }

  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  inline double Minimum(double a, double b)
  {
    return (a < b ? a: b);
  }

}

 // -----------------------------------------------------------------------------
 //
 // -----------------------------------------------------------------------------
 ReconstructionEngine::ReconstructionEngine()
 {

#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
    tbb::task_scheduler_init init;
    m_NumThreads = init.default_num_threads();
#else
    m_NumThreads = 1;
#endif
    setVerbose(true); //set this to enable cout::'s
    setVeryVerbose(true); //set this to ennable even more cout:: s
 }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionEngine::~ReconstructionEngine()
{
}


#include "ReconstructionEngine_Extra.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::InitializeTomoInputs(TomoInputsPtr v)
{
   v->sinoFile = "";
   v->initialReconFile = "";
   v->gainsInputFile = "";
   v->offsetsInputFile = "";
   v->varianceInputFile = "";
   v->InterpFlag=0;
   v->interpolateFactor=0.0;
   v->reconstructedOutputFile = "";
   v->tempDir = "";
   v->NumIter = 0;
   v->NumOuterIter = 0;
   v->SigmaX = 0.0;
   v->p = 0.0;
   v->StopThreshold = 0.0;
   v->useSubvolume = false;
   v->xStart = 0;
   v->xEnd = 0;
   v->yStart = 0;
   v->yEnd = 0;
   v->zStart = 0;
   v->zEnd = 0;
   v->fileXSize = 0;
   v->fileYSize = 0;
   v->fileZSize = 0;
   v->LengthZ = 0;
   v->delta_xz = 0;
   v->delta_xy = 0;
   v->tiltSelection = SOC::A_Tilt;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::InitializeSinogram(SinogramPtr v)
{
  v->N_r = 0;
  v->N_t = 0;
  v->N_theta = 0;
  v->delta_r = 0.0;
  v->delta_t = 0.0;
  v->counts = RealVolumeType::NullPointer();
  v->R0 = 0.0;
  v->RMax = 0.0;
  v->T0 = 0.0;
  v->TMax = 0.0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::InitializeGeometry(GeometryPtr v)
{
  v->Object = RealVolumeType::NullPointer();

  v->LengthX = 0.0;
  v->LengthY = 0.0;
  v->N_x = 0;
  v->N_z = 0;
  v->N_y = 0;
  v->x0 = 0.0;
  v->z0 = 0.0;
  v->y0 = 0.0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::InitializeAdvancedParams(AdvancedParametersPtr v)
{
  v->X_SHRINK_FACTOR = 0.6;
  v->X_STRETCH = 1;
  v->Z_STRETCH = 2;
  v->DETECTOR_RESPONSE_BINS = 64;
  v->PROFILE_RESOLUTION = 1536;
  v->BEAM_RESOLUTION = 512;
  v->AREA_WEIGHTED = 1;
  v->THRESHOLD_REDUCTION_FACTOR = 1;
  v->JOINT_ESTIMATION = 0;
  v->ZERO_SKIPPING = 1;
  v->NOISE_ESTIMATION = 1;

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::execute()
{
  uint64_t totalTime = EIMTOMO_getMilliSeconds();
  int32_t err = 0;

  std::stringstream ss;

  // int16_t i,j,k,Idx;
  size_t dims[3];

 // Int32ArrayType::Pointer Counter;
  UInt8Image_t::Pointer visitCount;

  uint16_t maxNumberOfDetectorElts;
  uint8_t status; //set to 1 if ICD has converged

  Real_t checksum = 0, temp;

  RealImageType::Pointer voxelProfile;
  RealVolumeType::Pointer detectorResponse;
  RealVolumeType::Pointer h_t;

  RealVolumeType::Pointer y_Est; //Estimated Sinogram
  RealVolumeType::Pointer finalSinogram; //To store and write the final sinogram resulting from our reconstruction
  RealVolumeType::Pointer errorSino; //Error Sinogram

  std::string indent("");

//#ifdef COST_CALCULATE //Commented out because if not the code fails to run.
  std::string filepath(m_TomoInputs->tempDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::CostFunctionFile);

  CostData::Pointer cost = CostData::New();
  cost->initOutputFile(filepath);
//#endif

#if OpenMBIR_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
#endif

  // Initialize the Sinogram
  if(m_TomoInputs == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The TomoInput Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }
  //Based on the inputs , calculate the "other" variables in the structure definition
  if(m_Sinogram == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The Sinogram Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }
  if(m_Geometry == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The Geometry Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }
  if (m_ForwardModel.get() == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The forward model was not set.", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }

  // Read the Input data from the supplied data file
  err = readInputData();
  if(err < 0)
  {
    return;
  }
  if (getCancel() == true) { setErrorCondition(-999); return; }

  err = m_ForwardModel->initializeBrightFieldData(m_Sinogram);
  if(err < 0)
  {
    return;
  }
  if (getCancel() == true) { setErrorCondition(-999); return; }


  err = m_ForwardModel->createNuisanceParameters(m_Sinogram);
  if(err < 0)
  {
    return;
  }

  if (getCancel() == true) { setErrorCondition(-999); return; }


  m_ForwardModel->printNuisanceParameters(m_Sinogram);

#ifdef BF_RECON //Take log of the input data after subtracting offset
    m_ForwardModel->processRawCounts(m_Sinogram);
#endif


  // Initialize the Geometry data from a rough reconstruction
  err = initializeRoughReconstructionData();
  if(err < 0)
  {
    return;
  }
  if (getCancel() == true) { setErrorCondition(-999); return; }

  // Get the actual boundaries in the X Direction since we "Extend Object" which makes
  // the output mrc file much wider than they really need to be.
  uint16_t cropStart=0;
  uint16_t cropEnd=m_Geometry->N_x;
  computeOriginalXDims(cropStart, cropEnd);
//  std::cout << "Crop Start: " << cropStart << std::endl;
//  std::cout << "Crop End:   " << cropEnd << std::endl;

  m_ForwardModel->allocateNuisanceParameters(m_Sinogram);

#ifdef COST_CALCULATE
  dims[0] = (m_TomoInputs->NumIter+1)*m_TomoInputs->NumOuterIter*3;
#endif

  dims[0] = m_Sinogram->N_theta;
  dims[1] = m_Sinogram->N_r;
  dims[2] = m_Sinogram->N_t;

  y_Est = RealVolumeType::New(dims, "y_Est");
  errorSino = RealVolumeType::New(dims, "ErrorSino");
  m_ForwardModel->weightInitialization(dims);
  finalSinogram = RealVolumeType::New(dims, "Final Sinogram");

  //calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram
  // Structure are converted to radians
  voxelProfile = calculateVoxelProfile(); //Verified with ML

  //Pre compute sine and cos theta to speed up computations
  HAADFDetectorParameters::Pointer haadfParameters = HAADFDetectorParameters::New();
  haadfParameters->setOffsetR ( ((m_TomoInputs->delta_xz / sqrt(3.0)) + m_Sinogram->delta_r / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS);
  haadfParameters->setOffsetT ( ((m_TomoInputs->delta_xz / 2) + m_Sinogram->delta_t / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS);
  haadfParameters->setBeamWidth(m_Sinogram->delta_r);
  haadfParameters->calculateSinCos(m_Sinogram);
  //Initialize the e-beam
  haadfParameters->initializeBeamProfile(m_Sinogram, m_AdvParams); //verified with ML


  // Initialize the Prior Model parameters - here we are using a QGGMRF Prior Model
  QGGMRF::QGGMRF_Values qggmrf_values;
  QGGMRF::initializePriorModel(m_TomoInputs, &qggmrf_values);
  // Set the prior model parameters into the Forward Model
  m_ForwardModel->setQGGMRFValues(&qggmrf_values);

  //globals assosiated with finding the optimal gain and offset parameters

  dims[0] = m_Sinogram->N_theta;
  dims[1] = 3;
  dims[2] = 0;
  //Hold the coefficients of a quadratic equation

  m_ForwardModel->costInitialization(m_Sinogram);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  //calculate sine and cosine of all angles and store in the global arrays sine and cosine
  DetectorResponse::Pointer dResponseFilter = DetectorResponse::New();
  dResponseFilter->setTomoInputs(m_TomoInputs);
  dResponseFilter->setSinogram(m_Sinogram);
  dResponseFilter->setAdvParams(m_AdvParams);
  dResponseFilter->setDetectorParameters(haadfParameters);
  dResponseFilter->setVoxelProfile(voxelProfile);
  dResponseFilter->setObservers(getObservers());
  dResponseFilter->setVerbose(getVerbose());
  dResponseFilter->setVeryVerbose(getVeryVerbose());
  dResponseFilter->execute();
  if(dResponseFilter->getErrorCondition() < 0)
  {
    std::cout << "Error Calling function detectorResponse in file " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
    setErrorCondition(-2);
    return;
  }
  detectorResponse = dResponseFilter->getResponse();
  // Writer the Detector Response to an output file
  DetectorResponseWriter::Pointer responseWriter = DetectorResponseWriter::New();
  responseWriter->setTomoInputs(m_TomoInputs);
  responseWriter->setSinogram(m_Sinogram);
  responseWriter->setAdvParams(m_AdvParams);
  responseWriter->setObservers(getObservers());
  responseWriter->setResponse(detectorResponse);
  responseWriter->setVerbose(getVerbose());
  responseWriter->setVeryVerbose(getVeryVerbose());
  responseWriter->execute();
  if(responseWriter->getErrorCondition() < 0)
  {
    std::cout << __FILE__ << "(" << __LINE__ << ") " << "Error writing detector response to file." <<  std::endl;
    setErrorCondition(-2);
    notify("Error Encountered During Reconstruction", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }

#ifdef RANDOM_ORDER_UPDATES
  dims[0] = m_Geometry->N_z;
  dims[1] = m_Geometry->N_x;
  dims[2] = 0;
  visitCount = UInt8Image_t::New(dims, "visitCount");
// Initialize the Array to zero
  ::memset(visitCount->d, 0, dims[0] * dims[1] * sizeof(uint8_t));
#endif//Random update


  //Gain and Offset Parameters Initialization
  m_ForwardModel->gainAndOffsetInitialization(m_Sinogram->N_theta);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  // Initialize H_t volume
  dims[0] = 1;
  dims[1] = m_Sinogram->N_theta;
  dims[2] = m_AdvParams->DETECTOR_RESPONSE_BINS;
  h_t = RealVolumeType::New(dims, "h_t");
  initializeHt(h_t, haadfParameters->getOffsetT() );

  checksum = 0;

  ss.str("");
  ss << "Calculating A Matrix....";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);



  std::vector<HAADFAMatrixCol::Pointer> voxelLineResponse(m_Geometry->N_y);

  maxNumberOfDetectorElts = (uint16_t)((m_TomoInputs->delta_xy / m_Sinogram->delta_t) + 2);
  dims[0] = maxNumberOfDetectorElts;
  for (uint16_t i = 0; i < m_Geometry->N_y; i++)
  {
    HAADFAMatrixCol::Pointer vlr = HAADFAMatrixCol::New(dims, 0);
    voxelLineResponse[i] = vlr;
  }

  //Calculating A-Matrix one column at a time
  //For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
  // And then store only the non zero entries by allocating a new array of the desired size
  std::vector<HAADFAMatrixCol::Pointer> tempCol(m_Geometry->N_x * m_Geometry->N_z);

  checksum = 0;
  temp = 0;
  uint32_t voxel_count = 0;
  for (uint16_t z = 0; z < m_Geometry->N_z; z++)
  {
    for (uint16_t x = 0; x < m_Geometry->N_x; x++)
    {
      tempCol[voxel_count] = HAADFAMatrixCol::calculateHAADFAMatrixColumnPartial(m_Sinogram, m_Geometry, m_TomoInputs, m_AdvParams,
                                                                                 z, x, 0, detectorResponse, haadfParameters);
      temp += tempCol[voxel_count]->count;
      if(0 == tempCol[voxel_count]->count )
      {
        //If this line is never hit and the Object is badly initialized
        //set it to zero
        for (uint16_t y = 0; y < m_Geometry->N_y; y++)
        {
          m_Geometry->Object->setValue(0.0, z, x, y);
        }
      }
      voxel_count++;
    }
  }

  storeVoxelResponse(h_t, voxelLineResponse, haadfParameters);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  if(getVerbose())
  {
    printf("Number of non zero entries of the forward projector is %lf\n", temp);
    printf("Geometry-Z %d\n", m_Geometry->N_z);
  }

  initializeVolume(y_Est, 0.0);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  //Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
  // Forward Project using the Forward Model
  err = m_ForwardModel->forwardProject(m_Sinogram, m_Geometry, tempCol, voxelLineResponse, y_Est , errorSino);
  if (err < 0)
  {
    return;
  }
  if (getCancel() == true) { setErrorCondition(-999); return; }


#ifdef FORWARD_PROJECT_MODE
  return 0; //exit the program once we finish forward projecting the object
#endif//Forward Project mode

#ifdef COST_CALCULATE
	err = calculateCost(cost,m_Sinogram,m_Geometry,errorSino,&qggmrf_values);
#endif //Cost calculation endif
	
//  int totalLoops = m_TomoInputs->NumOuterIter * m_TomoInputs->NumIter;
  //Loop through every voxel updating it by solving a cost function
  for (int16_t reconOuterIter = 0; reconOuterIter < m_TomoInputs->NumOuterIter; reconOuterIter++)
  {
    ss.str(""); // Clear the string stream
    indent = "";

    //The first time we may need to update voxels multiple times and then on just optimize over I,d,sigma,f once each outer loop
    if(reconOuterIter != 0)
    {
      m_TomoInputs->NumIter = 1;
    }

    for (int16_t reconInnerIter = 0; reconInnerIter < m_TomoInputs->NumIter; reconInnerIter++)
    {
      ss.str("");
      ss << "Outer Iterations: " << reconOuterIter << "/" << m_TomoInputs->NumOuterIter << " Inner Iterations: " << reconInnerIter << "/" << m_TomoInputs->NumIter << std::endl;

      indent = "    ";
      // This is all done PRIOR to calling what will become a method


      // This could contain multiple Subloops also
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      status = m_ForwardModel->updateVoxels(m_Sinogram, m_Geometry,
                                            reconOuterIter, reconInnerIter,
                                             visitCount, tempCol,
                                             errorSino, voxelLineResponse,
                                             cost);
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      if(status == 0)
      {
        break; //stop inner loop if we have hit the threshold value for x
      }
      // Check to see if we are canceled.
      if (getCancel() == true) { setErrorCondition(-999); return; }

      // Write out the VTK file
  //    {
  //      ss.str("");
  //      ss << m_TomoInputs->tempDir << MXADir::getSeparator() << reconOuterIter << "_" << ScaleOffsetCorrection::ReconstructedVtkFile;
  //      writeVtkFile(ss.str());
  //    }
      // Write out the MRC File
      {
        ss.str("");
        ss << m_TomoInputs->tempDir << MXADir::getSeparator() << reconOuterIter << "_" << reconInnerIter << "_" << ScaleOffsetCorrection::ReconstructedMrcFile;
        writeMRCFile(ss.str(), cropStart, cropEnd);
        notify(ss.str(), 0, Observable::UpdateIntermediateImage);
        m_TomoInputs->tempFiles.push_back(ss.str());
      }
#ifdef COST_CALCULATE
		
		/*********************Cost Calculation*************************************/
		int16_t err = calculateCost(cost,m_Sinogram,m_Geometry,errorSino,&qggmrf_values);
		if(err < 0)
		{
			std::cout<<"Cost went up after gain+offset update"<<std::endl;
			break;
		}
		/**************************************************************************/
#endif //Cost calculation endif  

    } /* ++++++++++ END Inner Iteration Loop +++++++++++++++ */
	  
	  
	if(0 == status && reconOuterIter >= 1) //
	{
		std::cout << "Exiting the code because status =0" << std::endl;
		break;
	}
	  
    if(m_AdvParams->JOINT_ESTIMATION)
    {
      m_ForwardModel->jointEstimation(m_Sinogram, errorSino, y_Est, cost);
#ifdef COST_CALCULATE
		//err = calculateCost(cost, Weight, errorSino);
		int16_t err = calculateCost(cost,m_Sinogram,m_Geometry,errorSino,&qggmrf_values);
		
		if(err < 0)
		{
			std::cout<<"Cost went up after gain+offset update"<<std::endl;
			break;
		}
#endif//cost
	
    } //Joint estimation endif

    if(m_AdvParams->NOISE_ESTIMATION)
    {
      m_ForwardModel->updateWeights(m_Sinogram, errorSino);
#ifdef COST_CALCULATE
		//err = calculateCost(cost, Weight, errorSino);
		int16_t err = calculateCost(cost,m_Sinogram,m_Geometry,errorSino,&qggmrf_values);
		if (err < 0)
		{
			std::cout<<"Cost went up after variance update"<<std::endl;
			break;
		}
#endif//cost
    }
    /*else
    {
      if(0 == status && reconOuterIter >= 1)
      { 
        std::cout << "Exiting the code because status =0" << std::endl;
        break;
      }
    } //Noise Model
    */

  }/* ++++++++++ END Outer Iteration Loop +++++++++++++++ */

  indent = "";
#if DEBUG_COSTS
  cost->printCosts(std::cout);
#endif

#ifdef FORWARD_PROJECT_MODE
  int fileError=0;
  FILE *Fp6;
  MAKE_OUTPUT_FILE(Fp6, fileError, m_TomoInputs->tempDir, ScaleOffsetCorrection::ForwardProjectedObjectFile);
  if (fileError == 1)
  {

  }
#endif

  if (getCancel() == true) { setErrorCondition(-999); return; }

  /* Write the Gains and Offsets to an output file */
  m_ForwardModel->writeNuisanceParameters(m_Sinogram);



  Real_t temp_final = 0.0;
  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
  {
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        temp_final = m_Sinogram->counts->getValue(i_theta, i_r, i_t) - errorSino->getValue(i_theta, i_r, i_t);
        finalSinogram->setValue(temp_final, i_theta, i_r, i_t);
      }
    }
  }

  if (getCancel() == true) { setErrorCondition(-999); return; }

 // This is writing the "ReconstructedSinogram.bin" file
  m_ForwardModel->writeSinogramFile(m_Sinogram, finalSinogram); // Writes the sinogram to a file

  // Writes ReconstructedObject.bin file
  {
      std::stringstream ss;
      ss << m_TomoInputs->tempDir << MXADir::getSeparator() << ScaleOffsetCorrection::ReconstructedObjectFile;
      writeReconstructionFile(ss.str());
  }
  // Write out the VTK file
  if (m_TomoInputs->vtkOutputFile.empty() == false)
  {
  //  std::stringstream ss;
  //  ss << m_TomoInputs->tempDir << MXADir::getSeparator() << ScaleOffsetCorrection::ReconstructedVtkFile;
    writeVtkFile(m_TomoInputs->vtkOutputFile, cropStart, cropEnd);
  }
  // Write out the MRC File
  if (m_TomoInputs->mrcOutputFile.empty() == false)
  {
  //  std::stringstream ss;
  //  ss << m_TomoInputs->tempDir << MXADir::getSeparator() << ScaleOffsetCorrection::ReconstructedMrcFile;
    writeMRCFile(m_TomoInputs->mrcOutputFile, cropStart, cropEnd);
  }

 // std::cout << "Should be writing .am file....  '" << m_TomoInputs->avizoOutputFile << "'"  << std::endl;
  if (m_TomoInputs->avizoOutputFile.empty() == false)
  {
    writeAvizoFile(m_TomoInputs->avizoOutputFile, cropStart, cropEnd);
  }

  std::cout << "Final Dimensions of Object: " << std::endl;
  std::cout << "  Nx = " << m_Geometry->N_x << std::endl;
  std::cout << "  Ny = " << m_Geometry->N_y << std::endl;
  std::cout << "  Nz = " << m_Geometry->N_z << std::endl;

  notify("Reconstruction Complete", 100, Observable::UpdateProgressValueAndMessage);
  setErrorCondition(0);
  std::cout << "Total Running Time for Execute: " << (EIMTOMO_getMilliSeconds() - totalTime) / 1000 << std::endl;
  return;
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RealImageType::Pointer ReconstructionEngine::calculateVoxelProfile()
{
  Real_t angle, maxValLineIntegral;
  Real_t temp, dist1, dist2, leftCorner, leftNear, rightNear, rightCorner, t;
  size_t dims[2] =
  { m_Sinogram->N_theta, m_AdvParams->PROFILE_RESOLUTION };
  RealImageType::Pointer voxProfile = RealImageType::New(dims, "VoxelProfile");

  Real_t checksum = 0;
  uint16_t i, j;
  FILE* fp = NULL;
  MAKE_OUTPUT_FILE(fp, m_TomoInputs->tempDir, ScaleOffsetCorrection::VoxelProfileFile);
  if(errno > 0)
  {
    std::string filepath(m_TomoInputs->tempDir);
    filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::VoxelProfileFile);\
    std::cout << "VoxelProfile will NOT be written to file '" << filepath << std::endl;
  }

  for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    m_Sinogram->angles[i] = m_Sinogram->angles[i] * (M_PI / 180.0);
    angle = m_Sinogram->angles[i];
    while (angle > M_PI_2)
      angle -= M_PI_2;

    while (angle < 0)
      angle += M_PI_2;

    if(angle <= M_PI_4)
    {
      maxValLineIntegral = m_TomoInputs->delta_xz / cos(angle);
    }
    else
    {
      maxValLineIntegral = m_TomoInputs->delta_xz / cos(M_PI_2 - angle);
    }
    temp = cos(M_PI_4);
    dist1 = temp * cos((M_PI_4 - angle));
    dist2 = temp * fabs((cos((M_PI_4 + angle))));
    leftCorner = 1 - dist1;
    leftNear = 1 - dist2;
    rightNear = 1 + dist2;
    rightCorner = 1 + dist1;

    for (j = 0; j < m_AdvParams->PROFILE_RESOLUTION; j++)
    {
      t = 2.0 * j / m_AdvParams->PROFILE_RESOLUTION; //2 is the normalized length of the profile (basically equl to 2*delta_xz)
      if(t <= leftCorner || t >= rightCorner) {voxProfile->setValue(0, i, j);}
      else if(t > rightNear) { voxProfile->setValue(maxValLineIntegral * (rightCorner - t) / (rightCorner - rightNear), i, j);}
      else if(t >= leftNear) { voxProfile->setValue(maxValLineIntegral, i, j);}
      else { voxProfile->setValue(maxValLineIntegral * (t - leftCorner) / (leftNear - leftCorner), i, j);}

      if(fp != NULL)
      {
        fwrite(voxProfile->getPointer(i, j), sizeof(Real_t), 1, fp);
      }
      checksum += voxProfile->getValue(i, j);
    }

  }

  //printf("Pixel Profile Check sum =%lf\n",checksum);
  if(fp != NULL)
  {
    fclose(fp);
  }
  return voxProfile;
}



// -----------------------------------------------------------------------------
//Finds the maximum of absolute value elements in an array
// -----------------------------------------------------------------------------
Real_t ReconstructionEngine::absMaxArray(std::vector<Real_t> &Array)
{
  uint16_t i;
  Real_t max;
  max = fabs(Array[0]);
  for(i =1; i < Array.size();i++)
    if(fabs(Array[i]) > max)
      max=fabs(Array[i]);
  return max;

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ReconstructionEngine::calculateCost(CostData::Pointer cost,
										SinogramPtr sinogram,
										GeometryPtr geometry,
										RealVolumeType::Pointer ErrorSino,
										QGGMRF::QGGMRF_Values* qggmrf_Values)
{
	Real_t cost_value = computeCost(sinogram, geometry, ErrorSino, qggmrf_Values);
	std::cout << "cost_value: " << cost_value << std::endl;
	int increase = cost->addCostValue(cost_value);
	if(increase == 1)
	{
		return -1;
	}
	cost->writeCostValue(cost_value);
	return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t ReconstructionEngine::computeCost(SinogramPtr sinogram, GeometryPtr geometry, RealVolumeType::Pointer ErrorSino, QGGMRF::QGGMRF_Values* qggmrf_values)
{
	Real_t cost = 0, temp = 0;
	Real_t delta;
	Real_t errSinoValue = 0.0;
	
	cost = m_ForwardModel->forwardCost(sinogram,ErrorSino);//Data term error
	cost += QGGMRF::PriorModelCost(geometry, qggmrf_values);//Prior model error

	return cost;
}


