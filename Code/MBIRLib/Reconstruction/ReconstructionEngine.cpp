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
#define PRINT_TIME(msg)\
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
   initVariables();
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
    tbb::task_scheduler_init init;
    m_NumThreads = init.default_num_threads();
#else
    m_NumThreads = 1;
#endif
    setVerbose(false); //set this to enable cout::'s
    setVeryVerbose(false); //set this to ennable even more cout:: s
 }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionEngine::~ReconstructionEngine()
{
}

 // These files are just Factored out CPP code because this file was getting really long
#include "ReconstructionEngine_UpdateVoxels.cpp"
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
  v->JOINT_ESTIMATION = 1;
  v->ZERO_SKIPPING = 1;
  v->NOISE_MODEL = 1;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::initVariables()
{
//  //Hamming Window here
//  HAMMING_WINDOW[0][0]= 0.0013; HAMMING_WINDOW[0][1]=0.0086; HAMMING_WINDOW[0][2]=0.0159; HAMMING_WINDOW[0][3]=0.0086;HAMMING_WINDOW[0][4]=0.0013;
//  HAMMING_WINDOW[1][0]= 0.0086; HAMMING_WINDOW[1][1]=0.0581;HAMMING_WINDOW[1][2]=0.1076;HAMMING_WINDOW[1][3]=0.0581;HAMMING_WINDOW[1][4]=0.0086;
//  HAMMING_WINDOW[2][0]= 0.0159;HAMMING_WINDOW[2][1]=0.1076;HAMMING_WINDOW[2][2]=0.1993;HAMMING_WINDOW[2][3]=0.1076;HAMMING_WINDOW[2][4]=0.0159;
//  HAMMING_WINDOW[3][0]= 0.0013;HAMMING_WINDOW[3][1]=0.0086;HAMMING_WINDOW[3][2]=0.0159;HAMMING_WINDOW[3][3]=0.0086;HAMMING_WINDOW[3][4]=0.0013;
//  HAMMING_WINDOW[4][0]= 0.0086;HAMMING_WINDOW[4][1]=0.0581;HAMMING_WINDOW[4][2]=0.1076;HAMMING_WINDOW[4][3]=0.0581;HAMMING_WINDOW[4][4]=0.0086;
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

  Int32ArrayType::Pointer Counter;
  UInt8Image_t::Pointer VisitCount;

  uint16_t MaxNumberOfDetectorElts;
  uint8_t status; //set to 1 if ICD has converged

  Real_t checksum = 0, temp;

  RealImageType::Pointer VoxelProfile;
  RealVolumeType::Pointer detectorResponse;
  RealVolumeType::Pointer H_t;

  RealVolumeType::Pointer Y_Est; //Estimated Sinogram
  RealVolumeType::Pointer Final_Sinogram; //To store and write the final sinogram resulting from our reconstruction
  RealVolumeType::Pointer ErrorSino; //Error Sinogram

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
    processRawCounts();
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

  Y_Est = RealVolumeType::New(dims, "Y_Est");
  ErrorSino = RealVolumeType::New(dims, "ErrorSino");
  m_ForwardModel->weightInitialization(dims);
  Final_Sinogram = RealVolumeType::New(dims, "Final Sinogram");

  //calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram
  // Structure are converted to radians
  VoxelProfile = calculateVoxelProfile(); //Verified with ML

  //Pre compute sine and cos theta to speed up computations
  HAADFDetectorParameters::Pointer haadfParameters = HAADFDetectorParameters::New();
  haadfParameters->setOffsetR ( ((m_TomoInputs->delta_xz / sqrt(3.0)) + m_Sinogram->delta_r / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS);
  haadfParameters->setOffsetT ( ((m_TomoInputs->delta_xz / 2) + m_Sinogram->delta_t / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS);
  haadfParameters->setBEAM_WIDTH(m_Sinogram->delta_r);
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
  NumOfViews = m_Sinogram->N_theta;
  LogGain = m_Sinogram->N_theta * log(m_ForwardModel->getTargetGain());

  m_ForwardModel->costInitialization(m_Sinogram);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  //calculate sine and cosine of all angles and store in the global arrays sine and cosine
  DetectorResponse::Pointer dResponseFilter = DetectorResponse::New();
  dResponseFilter->setTomoInputs(m_TomoInputs);
  dResponseFilter->setSinogram(m_Sinogram);
  dResponseFilter->setAdvParams(m_AdvParams);
  dResponseFilter->setDetectorParameters(haadfParameters);
  dResponseFilter->setVoxelProfile(VoxelProfile);
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
  VisitCount = UInt8Image_t::New(dims, "VisitCount");
// Initialize the Array to zero
  ::memset(VisitCount->d, 0, dims[0] * dims[1] * sizeof(uint8_t));
#endif//Random update


  //Gain and Offset Parameters Initialization
  m_ForwardModel->gainAndOffsetInitialization(m_Sinogram->N_theta);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  // Initialize H_t volume
  dims[0] = 1;
  dims[1] = m_Sinogram->N_theta;
  dims[2] = m_AdvParams->DETECTOR_RESPONSE_BINS;
  H_t = RealVolumeType::New(dims, "H_t");
  initializeHt(H_t, haadfParameters->getOffsetT() );

  checksum = 0;

  ss.str("");
  ss << "Calculating A Matrix....";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);



  std::vector<HAADFAMatrixCol::Pointer> VoxelLineResponse(m_Geometry->N_y);

  MaxNumberOfDetectorElts = (uint16_t)((m_TomoInputs->delta_xy / m_Sinogram->delta_t) + 2);
  dims[0] = MaxNumberOfDetectorElts;
  for (uint16_t i = 0; i < m_Geometry->N_y; i++)
  {
    HAADFAMatrixCol::Pointer vlr = HAADFAMatrixCol::New(dims, 0);
    VoxelLineResponse[i] = vlr;
  }

  //Calculating A-Matrix one column at a time
  //For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
  // And then store only the non zero entries by allocating a new array of the desired size
  std::vector<HAADFAMatrixCol::Pointer> TempCol(m_Geometry->N_x * m_Geometry->N_z);

  checksum = 0;
  temp = 0;
  uint32_t voxel_count = 0;
  for (uint16_t z = 0; z < m_Geometry->N_z; z++)
  {
    for (uint16_t x = 0; x < m_Geometry->N_x; x++)
    {
      TempCol[voxel_count] = HAADFAMatrixCol::calculateHAADFAMatrixColumnPartial(m_Sinogram, m_Geometry, m_TomoInputs, m_AdvParams,
                                                                                 z, x, 0, detectorResponse, haadfParameters);
      temp += TempCol[voxel_count]->count;
      if(0 == TempCol[voxel_count]->count )
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

  storeVoxelResponse(H_t, VoxelLineResponse, haadfParameters);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  if(getVerbose())
  {
    printf("Number of non zero entries of the forward projector is %lf\n", temp);
    printf("Geometry-Z %d\n", m_Geometry->N_z);
  }

  //Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
  //is Y_Est initailized to zero?
  initializeVolume(Y_Est, 0.0);

  if (getCancel() == true) { setErrorCondition(-999); return; }

  // Forward Project using the Forward Model
  err = m_ForwardModel->forwardProject(m_Sinogram, m_Geometry, TempCol, VoxelLineResponse, Y_Est , ErrorSino);
  if (err < 0)
  {
    return;
  }
  if (getCancel() == true) { setErrorCondition(-999); return; }


#ifdef FORWARD_PROJECT_MODE
  return 0; //exit the program once we finish forward projecting the object
#endif//Forward Project mode

#ifdef COST_CALCULATE
  err = calculateCost(cost, Weight, ErrorSino);
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
                                             VisitCount, TempCol,
                                             ErrorSino, VoxelLineResponse,
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

    } /* ++++++++++ END Inner Iteration Loop +++++++++++++++ */

    if(m_AdvParams->JOINT_ESTIMATION)
    {
      err = m_ForwardModel->jointEstimation(m_Sinogram, ErrorSino, Y_Est, cost);
      if(err < 0)
      {
        break;
      }
    } //Joint estimation endif

    if(m_AdvParams->NOISE_MODEL)
    {
      m_ForwardModel->updateWeights(m_Sinogram, ErrorSino);
#ifdef COST_CALCULATE
      err = calculateCost(cost, Weight, ErrorSino);
      if (err < 0)
      {
        std::cout<<"Cost went up after variance update"<<std::endl;
        break;
      }
#endif//cost
      if(0 == status && reconOuterIter >= 1) //&& VarRatio < STOPPING_THRESHOLD_Var_k && I_kRatio < STOPPING_THRESHOLD_I_k && Delta_kRatio < STOPPING_THRESHOLD_Delta_k)
      {
        std::cout << "Exiting the code because status =0" << std::endl;
        break;
      }
    }
    else
    {
      if(0 == status && reconOuterIter >= 1)
      { //&& I_kRatio < STOPPING_THRESHOLD_I_k && Delta_kRatio < STOPPING_THRESHOLD_Delta_k)
        std::cout << "Exiting the code because status =0" << std::endl;
        break;
      }
    } //Noise Model


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
        temp_final = m_Sinogram->counts->getValue(i_theta, i_r, i_t) - ErrorSino->getValue(i_theta, i_r, i_t);
        Final_Sinogram->setValue(temp_final, i_theta, i_r, i_t);
      }
    }
  }

  if (getCancel() == true) { setErrorCondition(-999); return; }

 // This is writing the "ReconstructedSinogram.bin" file
  m_ForwardModel->writeSinogramFile(m_Sinogram, Final_Sinogram); // Writes the sinogram to a file

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


/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void ReconstructionEngine::minMax(Real_t *low,Real_t *high, Real_t currentVoxelValue)
{
  uint8_t i,j,k;
  *low=NEIGHBORHOOD[INDEX_3(0,0,0)];
  *high=NEIGHBORHOOD[INDEX_3(0,0,0)];

  for(i = 0; i < 3;i++)
  {
    for(j=0; j < 3; j++)
    {
      for(k = 0; k < 3; k++)
      {
        //  if(NEIGHBORHOOD[i][j][k] != 0)
        //  printf("%lf ", NEIGHBORHOOD[i][j][k]);

        if(NEIGHBORHOOD[INDEX_3(i,j,k)] < *low)
          *low = NEIGHBORHOOD[INDEX_3(i,j,k)];
        if(NEIGHBORHOOD[INDEX_3(i,j,k)] > *high)
          *high=NEIGHBORHOOD[INDEX_3(i,j,k)];
      }
      //  printf("\n");
    }
  }


  if(THETA2 !=0)
  {
  *low = (*low > (currentVoxelValue - (THETA1/THETA2)) ? (currentVoxelValue - (THETA1/THETA2)): *low);

  *high = (*high < (currentVoxelValue - (THETA1/THETA2)) ? (currentVoxelValue - (THETA1/THETA2)): *high);
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RealImageType::Pointer ReconstructionEngine::calculateVoxelProfile()
{
    Real_t angle,MaxValLineIntegral;
    Real_t temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;
    size_t dims[2] = {m_Sinogram->N_theta ,m_AdvParams->PROFILE_RESOLUTION};
    RealImageType::Pointer VoxProfile = RealImageType::New(dims, "VoxelProfile");

    Real_t checksum=0;
    uint16_t i,j;
    FILE* Fp = NULL;
    MAKE_OUTPUT_FILE(Fp, m_TomoInputs->tempDir, ScaleOffsetCorrection::VoxelProfileFile);
    if (errno > 0)
    {
        std::string filepath(m_TomoInputs->tempDir);
        filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::VoxelProfileFile);\
        std::cout << "VoxelProfile will NOT be written to file '" << filepath << std::endl;
    }

    for (i=0;i<m_Sinogram->N_theta;i++)
    {
        m_Sinogram->angles[i]=m_Sinogram->angles[i]*(M_PI/180.0);
        angle=m_Sinogram->angles[i];
        while(angle > M_PI_2)
            angle -= M_PI_2;

        while(angle < 0)
            angle +=M_PI_2;

        if(angle <= M_PI_4)
        {
            MaxValLineIntegral = m_TomoInputs->delta_xz/cos(angle);
        }
        else
        {
            MaxValLineIntegral = m_TomoInputs->delta_xz/cos(M_PI_2-angle);
        }
        temp=cos(M_PI_4);
        dist1 = temp * cos((M_PI_4 - angle));
        dist2 = temp * fabs((cos((M_PI_4 + angle))));
        LeftCorner = 1-dist1;
        LeftNear = 1-dist2;
        RightNear = 1+dist2;
        RightCorner = 1+dist1;

        for(j = 0;j<m_AdvParams->PROFILE_RESOLUTION;j++)
        {
            t = 2.0*j /m_AdvParams->PROFILE_RESOLUTION;//2 is the normalized length of the profile (basically equl to 2*delta_xz)
            if(t <= LeftCorner || t >= RightCorner)
                VoxProfile->setValue(0, i, j);
            else if(t > RightNear)
                VoxProfile->setValue(MaxValLineIntegral*(RightCorner-t)/(RightCorner-RightNear), i, j);
            else if(t >= LeftNear)
                VoxProfile->setValue(MaxValLineIntegral, i, j);
            else
                VoxProfile->setValue(MaxValLineIntegral*(t-LeftCorner)/(LeftNear-LeftCorner), i, j);

            if (Fp != NULL)
            {
                fwrite( VoxProfile->getPointer(i, j), sizeof(Real_t),1,Fp);
            }
            checksum+=VoxProfile->getValue(i, j);
        }

    }

    //printf("Pixel Profile Check sum =%lf\n",checksum);
    if (Fp != NULL) {
        fclose(Fp);
    }
    return VoxProfile;
}

#if 0
/*******************************************************************
 Forwards Projects the Object and stores it in a 3-D matrix
 ********************************************************************/
RealVolumeType::Pointer ReconstructionEngine::forwardProject(RealVolumeType::Pointer DetectorResponse,
                                                  RealVolumeType::Pointer H_t)
{
  notify("Executing Forward Projection", 50, Observable::UpdateProgressValueAndMessage);

  Real_t x, z, y;
  Real_t r, rmin, rmax, t, tmin, tmax;
  Real_t w1, w2, f1, f2;
  Real_t center_r, delta_r, center_t, delta_t;
  int16_t index_min, index_max, slice_index_min, slice_index_max, i_r, i_t;

  int16_t index_delta_t, index_delta_r;
  size_t dims[3] =
  { m_Sinogram->N_theta, m_Sinogram->N_r, m_Sinogram->N_t };
  RealVolumeType::Pointer Y_Est = RealVolumeType::New(dims, "Y_Est");

  for (int16_t j = 0; j < m_Geometry->N_z; j++)
  {
    for (int16_t k = 0; k < m_Geometry->N_x; k++)
    {
      x = m_Geometry->x0 + ((Real_t)k + 0.5) * m_TomoInputs->delta_xz; //0.5 is for center of voxel. x_0 is the left corner
      z = m_Geometry->z0 + ((Real_t)j + 0.5) * m_TomoInputs->delta_xz; //0.5 is for center of voxel. z_0 is the left corner

      for (int16_t i = 0; i < m_Geometry->N_y; i++)
      {
        for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
        {
          r = x * cosine->d[i_theta] - z * sine->d[i_theta];
          rmin = r - m_TomoInputs->delta_xz;
          rmax = r + m_TomoInputs->delta_xz;

          if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax) continue;

          index_min = static_cast<uint16_t>(floor(((rmin - m_Sinogram->R0) / m_Sinogram->delta_r)));
          index_max = static_cast<uint16_t>(floor((rmax - m_Sinogram->R0) / m_Sinogram->delta_r));

          if(index_max >= m_Sinogram->N_r) index_max = m_Sinogram->N_r - 1;

          if(index_min < 0) index_min = 0;

          y = m_Geometry->y0 + ((double)i + 0.5) * m_TomoInputs->delta_xy;
          t = y;

          tmin = (t - m_TomoInputs->delta_xy / 2) > m_Sinogram->T0 ? t - m_TomoInputs->delta_xy / 2 : m_Sinogram->T0;
          tmax = (t + m_TomoInputs->delta_xy / 2) <= m_Sinogram->TMax ? t + m_TomoInputs->delta_xy / 2 : m_Sinogram->TMax;

          slice_index_min = static_cast<uint16_t>(floor((tmin - m_Sinogram->T0) / m_Sinogram->delta_t));
          slice_index_max = static_cast<uint16_t>(floor((tmax - m_Sinogram->T0) / m_Sinogram->delta_t));

          if(slice_index_min < 0) slice_index_min = 0;
          if(slice_index_max >= m_Sinogram->N_t) slice_index_max = m_Sinogram->N_t - 1;

          for (i_r = index_min; i_r <= index_max; i_r++)
          {
            center_r = m_Sinogram->R0 + ((double)i_r + 0.5) * m_Sinogram->delta_r;
            delta_r = fabs(center_r - r);
            index_delta_r = static_cast<uint16_t>(floor((delta_r / OffsetR)));

            if(index_delta_r < m_AdvParams->DETECTOR_RESPONSE_BINS)
            {
              w1 = delta_r - index_delta_r * OffsetR;
              w2 = (index_delta_r + 1) * OffsetR - delta_r;
              uint16_t iidx = index_delta_r + 1 < m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_r + 1 : m_AdvParams->DETECTOR_RESPONSE_BINS - 1;
              f1 = (w2 / OffsetR) * DetectorResponse->getValue(0, i_theta, index_delta_r)
                  +(w1 / OffsetR) * DetectorResponse->getValue(0, i_theta, iidx);
            }
            else
            {
              f1 = 0;
            }

            for (i_t = slice_index_min; i_t <= slice_index_max; i_t++)
            {
              center_t = m_Sinogram->T0 + ((double)i_t + 0.5) * m_Sinogram->delta_t;
              delta_t = fabs(center_t - t);
              index_delta_t = static_cast<uint16_t>(floor((delta_t / OffsetT)));

              if(index_delta_t < m_AdvParams->DETECTOR_RESPONSE_BINS)
              {
                w1 = delta_t - index_delta_t * OffsetT;
                w2 = (index_delta_t + 1) * OffsetT - delta_t;
                uint16_t iidx = index_delta_t + 1 < m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : m_AdvParams->DETECTOR_RESPONSE_BINS - 1;
                f2 = (w2 / OffsetT) * H_t->getValue(0, i_theta, index_delta_t) + (w1 / OffsetT) * H_t->getValue(0, i_theta, iidx);
              }
              else
              {
                f2 = 0;
              }
              //Y_Est->d[i_theta][i_r][i_t] += f1*f2*m_Geometry->Object->getValue(j, k, i);
              Y_Est->addToValue(f1 * f2 * m_Geometry->Object->getValue(j, k, i), i_theta, i_r, i_t);
            }
          }
        }
      }

    }

  }

  return Y_Est;
}
#endif


#ifndef EIMTOMO_USE_QGGMRF
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double ReconstructionEngine::surrogateFunctionBasedMin(Real_t currentVoxelValue)
{
  double numerator_sum = 0;
  double denominator_sum = 0;
  double alpha, update = 0;
  double product = 1;
  uint8_t i, j, k;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        if(currentVoxelValue != NEIGHBORHOOD[i][j][k])
        {
          product = ((double)FILTER[i][j][k] * pow(fabs(currentVoxelValue - NEIGHBORHOOD[i][j][k]), (MRF_P - 2.0)));
          numerator_sum += (product * (currentVoxelValue - NEIGHBORHOOD[i][j][k]));
          denominator_sum += product;
        }
      }
    }
  }
  numerator_sum /= SIGMA_X_P;
  denominator_sum /= SIGMA_X_P;

  numerator_sum += THETA1;
  denominator_sum += THETA2;

  if(THETA2 > 0)
  {
    alpha = (-1 * numerator_sum) / (denominator_sum);
    update = currentVoxelValue + Detail::Clip(alpha, -currentVoxelValue, std::numeric_limits<float>::infinity());
  }
  else
  {
    update = 0;
  }

  if(update > 70000) printf("%lf\n", update);

  return update;

}
#endif

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



