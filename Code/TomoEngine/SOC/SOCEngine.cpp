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
#include "SOCEngine.h"

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
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/Common/EIMTime.h"
#include "TomoEngine/Common/CE_ConstraintEquation.hpp"
#include "TomoEngine/Common/DerivOfCostFunc.hpp"
#include "TomoEngine/IO/RawGeometryWriter.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"
#include "TomoEngine/IO/MRCWriter.h"
#include "TomoEngine/IO/RawGeometryWriter.h"
#include "TomoEngine/IO/NuisanceParamWriter.h"
#include "TomoEngine/IO/NuisanceParamReader.h"
#include "TomoEngine/IO/SinogramBinWriter.h"
#include "TomoEngine/IO/VTKFileWriters.hpp"
#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/SOC/ForwardProject.h"

#include "TomoEngine/Filters/ComputeInitialOffsets.h"
#include "TomoEngine/Filters/DetectorResponse.h"
#include "TomoEngine/Filters/DetectorResponseWriter.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/Filters/MRCSinogramInitializer.h"
#include "TomoEngine/Filters/RawSinogramInitializer.h"
#include "TomoEngine/Filters/GainsOffsetsReader.h"
#include "TomoEngine/Filters/ComputeInitialOffsets.h"
#include "TomoEngine/Filters/InitialReconstructionInitializer.h"
#include "TomoEngine/Filters/InitialReconstructionBinReader.h"



#define USE_TBB_TASK_GROUP 1
#if defined (TomoEngine_USE_PARALLEL_ALGORITHMS)
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#include <tbb/task.h>

#endif

#define START_TIMER startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME(msg)\
    std::cout << indent << msg << ": " << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

#define MAKE_OUTPUT_FILE(Fp, err, outdir, filename)\
    {\
    std::string filepath(outdir);\
    filepath = filepath.append(MXADir::getSeparator()).append(filename);\
    errno = 0;\
    err = 0;\
    Fp = fopen(filepath.c_str(),"wb");\
    if (Fp == NULL) { std::cout << "Error " << errno << " Opening Output file " << filepath << std::endl; err = -1; }\
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
  //
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
 SOCEngine::SOCEngine()
 {
   initVariables();
#if defined (TomoEngine_USE_PARALLEL_ALGORITHMS)
    tbb::task_scheduler_init init;
    m_NumThreads = init.default_num_threads();
#else
    m_NumThreads = 1;
#endif
    //setVerbose(true); //set this to enable cout::'s
    //setVeryVerbose(true); //set this to ennable even more cout:: s
 }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SOCEngine::~SOCEngine()
{
}

 // These files are just Factored out CPP code because this file was getting really long
#include "SOCEngine_UpdateVoxels.cpp"
#include "SOCEngine_Extra.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::InitializeTomoInputs(TomoInputsPtr v)
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
   v->defaultOffset = 0.0;
   v->useDefaultOffset = false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::InitializeSinogram(SinogramPtr v)
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
  v->targetGain = 0.0;
  v->InitialGain = RealArrayType::NullPointer();
  v->InitialOffset = RealArrayType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::InitializeGeometry(GeometryPtr v)
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
void SOCEngine::InitializeScaleOffsetParams(ScaleOffsetParamsPtr v)
{
  v->I_0 = RealArrayType::NullPointer();
  v->mu = RealArrayType::NullPointer();
  v->alpha = RealArrayType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::InitializeAdvancedParams(AdvancedParametersPtr v)
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
void SOCEngine::initVariables()
{
  FILTER[INDEX_3(0,0,0)] = 0.0302; FILTER[INDEX_3(0,0,1)] = 0.0370; FILTER[INDEX_3(0,0,2)] = 0.0302;
  FILTER[INDEX_3(0,1,0)] = 0.0370; FILTER[INDEX_3(0,1,1)] = 0.0523; FILTER[INDEX_3(0,1,2)] = 0.0370;
  FILTER[INDEX_3(0,2,0)] = 0.0302; FILTER[INDEX_3(0,2,1)] = 0.0370; FILTER[INDEX_3(0,2,2)] = 0.0302;

  FILTER[INDEX_3(1,0,0)] = 0.0370; FILTER[INDEX_3(1,0,1)] = 0.0523; FILTER[INDEX_3(1,0,2)] = 0.0370;
  FILTER[INDEX_3(1,1,0)] = 0.0523; FILTER[INDEX_3(1,1,1)] = 0.0000; FILTER[INDEX_3(1,1,2)] = 0.0523;
  FILTER[INDEX_3(1,2,0)] = 0.0370; FILTER[INDEX_3(1,2,1)] = 0.0523; FILTER[INDEX_3(1,2,2)] = 0.0370;

  FILTER[INDEX_3(2,0,0)] = 0.0302; FILTER[INDEX_3(2,0,1)] = 0.0370; FILTER[INDEX_3(2,0,2)] = 0.0302;
  FILTER[INDEX_3(2,1,0)] = 0.0370; FILTER[INDEX_3(2,1,1)] = 0.0523; FILTER[INDEX_3(2,1,2)] = 0.0370;
  FILTER[INDEX_3(2,2,0)] = 0.0302; FILTER[INDEX_3(2,2,1)] = 0.0370; FILTER[INDEX_3(2,2,2)] = 0.0302;

  //Hamming Window here
  HAMMING_WINDOW[0][0]= 0.0013; HAMMING_WINDOW[0][1]=0.0086; HAMMING_WINDOW[0][2]=0.0159; HAMMING_WINDOW[0][3]=0.0086;HAMMING_WINDOW[0][4]=0.0013;
  HAMMING_WINDOW[1][0]= 0.0086; HAMMING_WINDOW[1][1]=0.0581;HAMMING_WINDOW[1][2]=0.1076;HAMMING_WINDOW[1][3]=0.0581;HAMMING_WINDOW[1][4]=0.0086;
  HAMMING_WINDOW[2][0]= 0.0159;HAMMING_WINDOW[2][1]=0.1076;HAMMING_WINDOW[2][2]=0.1993;HAMMING_WINDOW[2][3]=0.1076;HAMMING_WINDOW[2][4]=0.0159;
  HAMMING_WINDOW[3][0]= 0.0013;HAMMING_WINDOW[3][1]=0.0086;HAMMING_WINDOW[3][2]=0.0159;HAMMING_WINDOW[3][3]=0.0086;HAMMING_WINDOW[3][4]=0.0013;
  HAMMING_WINDOW[4][0]= 0.0086;HAMMING_WINDOW[4][1]=0.0581;HAMMING_WINDOW[4][2]=0.1076;HAMMING_WINDOW[4][3]=0.0581;HAMMING_WINDOW[4][4]=0.0086;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::execute()
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
  RealVolumeType::Pointer Weight; //This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation

  std::string indent("");

//#ifdef COST_CALCULATE //Commented out because if not the code fails to run.
  std::string filepath(m_TomoInputs->tempDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::CostFunctionFile);

  CostData::Pointer cost = CostData::New();
  cost->initOutputFile(filepath);
//#endif

#if TomoEngine_USE_PARALLEL_ALGORITHMS
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

  // Read the Input data from the supplied data file
  err = readInputData();
  if(err < 0)
  {
    return;
  }

  err = initializeBrightFieldData();
  if(err < 0)
  {
    return;
  }

  err = createInitialGainsData();
  if(err < 0)
  {
    return;
  }
  err = createInitialOffsetsData();
  if(err < 0)
  {
    return;
  }
  err = createInitialVariancesData();
  if(err < 0)
  {
    return;
  }

  if(getVeryVerbose())
  {
    // Print out the Initial Gains, Offsets, Variances
    std::cout << "---------------- Initial Gains, Offsets, Variances -------------------" << std::endl;
    std::cout << "Tilt\tGain\tOffset";

    if(NULL != m_Sinogram->InitialVariance.get())
    {
      std::cout << "\tVariance";
    }
    std::cout << std::endl;

    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      std::cout << i_theta << "\t" << m_Sinogram->InitialGain->d[i_theta] << "\t" << m_Sinogram->InitialOffset->d[i_theta];
      if(NULL != m_Sinogram->InitialVariance.get())
      {
        std::cout << "\t" << m_Sinogram->InitialVariance->d[i_theta];
      }
      std::cout << std::endl;
    }
  }

#ifdef BF_RECON //Take log of the input data after subtracting offset
	processRawCounts();
#endif
	
  // Initialize the Geometry data from a rough reconstruction
  err = initializeRoughReconstructionData();
  if(err < 0)
  {
    return;
  }


  //Gain, Offset and Variance Parameter Structures
  ScaleOffsetParamsPtr NuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);
  dims[1] = m_Sinogram->N_t;
  dims[0] = m_Sinogram->N_theta;
  NuisanceParams->I_0 = RealArrayType::New(dims, "NuisanceParams->I_0");
  NuisanceParams->mu = RealArrayType::New(dims, "NuisanceParams->mu");
  if(m_AdvParams->NOISE_MODEL)
  {
    //alpha is the noise variance adjustment factor
    NuisanceParams->alpha = RealArrayType::New(dims, "NuisanceParams->alpha");
  }
  else
  {
    NuisanceParams->alpha = RealArrayType::NullPointer();
  }

#if ROI
  UInt8Image_t::Pointer Mask;
//  DATA_TYPE EllipseA,EllipseB;
#endif

#ifdef COST_CALCULATE
  dims[0] = (m_TomoInputs->NumIter+1)*m_TomoInputs->NumOuterIter*3;
#endif

  dims[0] = m_Sinogram->N_theta;
  dims[1] = m_Sinogram->N_r;
  dims[2] = m_Sinogram->N_t;

  Y_Est = RealVolumeType::New(dims, "Y_Est");
  ErrorSino = RealVolumeType::New(dims, "ErrorSino");
  Weight = RealVolumeType::New(dims, "Weight");
  Final_Sinogram = RealVolumeType::New(dims, "Final Sinogram");

  //Setting the value of all the private members
  OffsetR = ((m_TomoInputs->delta_xz / sqrt(3.0)) + m_Sinogram->delta_r / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS;
  OffsetT = ((m_TomoInputs->delta_xz / 2) + m_Sinogram->delta_t / 2) / m_AdvParams->DETECTOR_RESPONSE_BINS;

  BEAM_WIDTH = m_Sinogram->delta_r;

#ifdef EIMTOMO_USE_QGGMRF
  m_QGGMRF_Values.MRF_P = 2;
  m_QGGMRF_Values.MRF_Q = m_TomoInputs->p;
  m_QGGMRF_Values.MRF_C = 0.01;
  m_QGGMRF_Values.MRF_ALPHA = 1.5;
  m_QGGMRF_Values.SIGMA_X_P = pow(m_TomoInputs->SigmaX, m_QGGMRF_Values.MRF_P);
  m_QGGMRF_Values.SIGMA_X_P_Q = pow(m_TomoInputs->SigmaX, (m_QGGMRF_Values.MRF_P - m_QGGMRF_Values.MRF_Q));
  m_QGGMRF_Values.SIGMA_X_Q = pow(m_TomoInputs->SigmaX, m_QGGMRF_Values.MRF_Q);
#else
  MRF_P = m_TomoInputs->p;
  SIGMA_X_P = pow(m_TomoInputs->SigmaX,MRF_P);
#endif //QGGMRF
  //globals assosiated with finding the optimal gain and offset parameters

  dims[0] = m_Sinogram->N_theta;
  dims[1] = 3;
  dims[2] = 0;
  //Hold the coefficients of a quadratic equation
  QuadraticParameters = RealImageType::New(dims, "QuadraticParameters");
  Qk_cost = RealImageType::New(dims, "Qk_cost");
  dims[1] = 2;
  bk_cost = RealImageType::New(dims, "bk_cost");

  dims[0] = m_Sinogram->N_theta;
  ck_cost = RealArrayType::New(dims, "ck_cost");

  NumOfViews = m_Sinogram->N_theta;
  LogGain = m_Sinogram->N_theta * log(m_Sinogram->targetGain);
  dims[0] = m_Sinogram->N_theta;
  d1 = RealArrayType::New(dims, "d1");
  d2 = RealArrayType::New(dims, "d2");

  //calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram Structure are converted to radians
  VoxelProfile = calculateVoxelProfile(); //Verified with ML
  //Pre compute sine and cos theta to speed up computations
  calculateSinCos();
  //Initialize the e-beam
  initializeBeamProfile(); //verified with ML

  //calculate sine and cosine of all angles and store in the global arrays sine and cosine
  DetectorResponse::Pointer dResponseFilter = DetectorResponse::New();
  dResponseFilter->setTomoInputs(m_TomoInputs);
  dResponseFilter->setSinogram(m_Sinogram);
  dResponseFilter->setAdvParams(m_AdvParams);
  dResponseFilter->setBeamWidth(BEAM_WIDTH);
  dResponseFilter->setOffsetR(OffsetR);
  dResponseFilter->setOffsetT(OffsetT);
  dResponseFilter->setVoxelProfile(VoxelProfile);
  dResponseFilter->setBeamProfile(BeamProfile);
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
    std::cout << "Error writing detector response to file." << __FILE__ << "(" << __LINE__ << ")" << std::endl;
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

  dims[0] = m_Geometry->N_z; //height
  dims[1] = m_Geometry->N_x; //width
  dims[2] = 0;

  MagUpdateMap = RealImageType::New(dims, "Update Map for voxel lines");
  FiltMagUpdateMap = RealImageType::New(dims, "Update Map for voxel lines");
  MagUpdateMask = UInt8Image_t::New(dims, "Update Mask for selecting voxel lines NHICD");

#if ROI
  //Mask = (uint8_t**)get_img(m_Geometry->N_x, m_Geometry->N_z,sizeof(uint8_t));//width,height
  dims[0] = m_Geometry->N_z;
  dims[1] = m_Geometry->N_x;
  Mask = UInt8Image_t::New(dims, "Mask");
  initializeROIMask(Mask);
#endif
  //m_Sinogram->targetGain=20000;

  //Gain and Offset Parameters Initialization
  gainAndOffsetInitialization(NuisanceParams);

  // Initialize H_t volume
  dims[0] = 1;
  dims[1] = m_Sinogram->N_theta;
  dims[2] = m_AdvParams->DETECTOR_RESPONSE_BINS;
  H_t = RealVolumeType::New(dims, "H_t");
  initializeHt(H_t);

  checksum = 0;

  ss.str("");
  ss << "Calculating AMatrix....";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);

  //TODO: All this needs to be deallocated at some point


  //AMatrixCol* VoxelLineResponse = (AMatrixCol*)get_spc(m_Geometry->N_y, sizeof(AMatrixCol));
  std::vector<AMatrixCol::Pointer> VoxelLineResponse(m_Geometry->N_y);

  MaxNumberOfDetectorElts = (uint16_t)((m_TomoInputs->delta_xy / m_Sinogram->delta_t) + 2);
  dims[0] = MaxNumberOfDetectorElts;
  for (uint16_t i = 0; i < m_Geometry->N_y; i++)
  {
    AMatrixCol::Pointer vlr = AMatrixCol::New(dims, 0);
    VoxelLineResponse[i] = vlr;
  }

  //Calculating A-Matrix one column at a time
  //For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
  // And then store only the non zero entries by allocating a new array of the desired size
  //AMatrixCol** TempCol = (AMatrixCol**)get_spc(m_Geometry->N_x * m_Geometry->N_z, sizeof(AMatrixCol*));
  std::vector<AMatrixCol::Pointer> TempCol(m_Geometry->N_x * m_Geometry->N_z);

  checksum = 0;
  temp = 0;
  uint32_t voxel_count = 0;
  for (uint16_t z = 0; z < m_Geometry->N_z; z++)
  {
    for (uint16_t x = 0; x < m_Geometry->N_x; x++)
    {
      TempCol[voxel_count] = calculateAMatrixColumnPartial(z, x, 0, detectorResponse);
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

  storeVoxelResponse(H_t, VoxelLineResponse);

  if(getVerbose())
  {
    printf("Number of non zero entries of the forward projector is %lf\n", temp);
    printf("Geometry-Z %d\n", m_Geometry->N_z);
  }

  //Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
  //is Y_Est initailized to zero?
  initializeVolume(Y_Est, 0.0);

  notify("Starting Forward Projection", 10, Observable::UpdateProgressValueAndMessage);
  START_TIMER;
  // This next section looks crazy with all the #if's but this makes sure we are
  // running the exact same code whether in parallel or serial.

#if TomoEngine_USE_PARALLEL_ALGORITHMS
  tbb::task_group* g = new tbb::task_group;
  if (getVerbose())
  {
    std::cout << "Default Number of Threads to Use: " << init.default_num_threads() << std::endl;
    std::cout << "Forward Projection Running in Parallel." << std::endl;
  }
#else
  if(getVerbose())
  {
    std::cout << "Forward Projection Running in Serial." << std::endl;
  }
#endif
  // Queue up a thread for each z layer of the Geometry. The threads will only be
  // run as hardware resources open up so this will not just fire up a gazillion
  // threads.
  for (uint16_t t = 0; t < m_Geometry->N_z; t++)
  {
#if TomoEngine_USE_PARALLEL_ALGORITHMS
    g->run(ForwardProject(m_Sinogram.get(), m_Geometry.get(), TempCol, VoxelLineResponse, Y_Est, NuisanceParams.get(), t, this));
#else
    ForwardProject fp(m_Sinogram.get(), m_Geometry.get(), TempCol, VoxelLineResponse, Y_Est, NuisanceParams.get(), t, this);
    fp.setObservers(getObservers());
    fp();
#endif
  }
#if TomoEngine_USE_PARALLEL_ALGORITHMS
  g->wait(); // Wait for all the threads to complete before moving on.
  delete g;
#endif

  STOP_TIMER;PRINT_TIME("Forward Project Time");

  //Calculate Error Sinogram - Can this be combined with previous loop?
  //Also compute weights of the diagonal covariance matrix
  calculateMeasurementWeight(Weight, NuisanceParams, ErrorSino, Y_Est);

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

    //The first time we may need to update voxels multiple times and then on just optimize over I,d,\sigma,f once each outer loop
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
      VoxelUpdateType updateType = RegularRandomOrderUpdate;
#ifdef NHICD
      if(0 == reconInnerIter % 2)
      {
        updateType = HomogeniousUpdate;
      }
      else
      {
        updateType = NonHomogeniousUpdate;
      }
#else

#endif//NHICD end if
      // This could contain multiple Subloops also
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      status =
          updateVoxels(reconOuterIter, reconInnerIter, updateType, VisitCount, TempCol, ErrorSino, Weight, VoxelLineResponse, NuisanceParams.get(), Mask, cost);
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      if(status == 0)
      {
        break; //stop inner loop if we have hit the threshold value for x
      }
      // Check to see if we are canceled.
      if(getCancel() == true)
      {
        setErrorCondition(-100);
        return;
      }

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
        writeMRCFile(ss.str());
        notify(ss.str(), 0, Observable::UpdateIntermediateImage);
      }

    } /* ++++++++++ END Inner Iteration Loop +++++++++++++++ */

	

    if(m_AdvParams->JOINT_ESTIMATION)
    {
      err = jointEstimation(Weight, NuisanceParams, ErrorSino, Y_Est, cost);
      if(err < 0)
      {
        break;
      }
    } //Joint estimation endif

	if(m_AdvParams->NOISE_MODEL)
	{
		  updateWeights(Weight, NuisanceParams, ErrorSino);
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
		  {//&& I_kRatio < STOPPING_THRESHOLD_I_k && Delta_kRatio < STOPPING_THRESHOLD_Delta_k)
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

  /* Write the Gains and Offsets to an output file */
  writeNuisanceParameters(NuisanceParams);

  if(getVerbose()) {
    std::cout << "Tilt\tFinal Gains\tFinal Offsets\tFinal Variances" << std::endl;
    for (uint16_t i_theta = 0; i_theta < getSinogram()->N_theta; i_theta++)
    {
	
		if(m_AdvParams->NOISE_MODEL)
		{
      std::cout << i_theta << "\t" << NuisanceParams->I_0->d[i_theta] <<
      "\t" << NuisanceParams->mu->d[i_theta] <<
      "\t" << NuisanceParams->alpha->d[i_theta] << std::endl;
		}
		else
		{
			std::cout << i_theta << "\t" << NuisanceParams->I_0->d[i_theta] <<
			"\t" << NuisanceParams->mu->d[i_theta] << std::endl;		
		}
    }
  }

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

  writeSinogramFile(NuisanceParams, Final_Sinogram); // Writes the sinogram to a file
  writeReconstructionFile(); // Writes the m_Geometry to a file
  // Write out the VTK file
  {
    std::stringstream ss;
    ss << m_TomoInputs->tempDir << MXADir::getSeparator() << ScaleOffsetCorrection::ReconstructedVtkFile;
    writeVtkFile(ss.str());
  }
  // Write out the MRC File
  {
    std::stringstream ss;
    ss << m_TomoInputs->tempDir << MXADir::getSeparator() << ScaleOffsetCorrection::ReconstructedMrcFile;
    writeMRCFile(ss.str());
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
void SOCEngine::minMax(Real_t *low,Real_t *high, Real_t currentVoxelValue)
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
RealImageType::Pointer SOCEngine::calculateVoxelProfile()
{
  Real_t angle,MaxValLineIntegral;
  Real_t temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;
  size_t dims[2] = {m_Sinogram->N_theta ,m_AdvParams->PROFILE_RESOLUTION};
  RealImageType::Pointer VoxProfile = RealImageType::New(dims, "VoxelProfile");

  Real_t checksum=0;
  uint16_t i,j;
  FILE* Fp = NULL;
  int fileError = 0;
  MAKE_OUTPUT_FILE(Fp, fileError, m_TomoInputs->tempDir, ScaleOffsetCorrection::VoxelProfileFile);
  if (fileError < 0)
  {

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

      fwrite( VoxProfile->getPointer(i, j), sizeof(Real_t),1,Fp);
      checksum+=VoxProfile->getValue(i, j);
    }

  }

  //printf("Pixel Profile Check sum =%lf\n",checksum);
  fclose(Fp);
  return VoxProfile;
}

/*******************************************************************
 Forwards Projects the Object and stores it in a 3-D matrix
 ********************************************************************/
RealVolumeType::Pointer SOCEngine::forwardProject(RealVolumeType::Pointer DetectorResponse,
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


/* Initializes the global variables cosine and sine to speed up computation
 */
void SOCEngine::calculateSinCos()
{
  uint16_t i;
  size_t dims[1] = { m_Sinogram->N_theta };
  cosine = RealArrayType::New(dims, "cosine");
  sine = RealArrayType::New(dims, "sine");

  for(i=0;i<m_Sinogram->N_theta;i++)
  {
    cosine->d[i]=cos(m_Sinogram->angles[i]);
    sine->d[i]=sin(m_Sinogram->angles[i]);
  }
}

void SOCEngine::initializeBeamProfile()
{
  uint16_t i;
  Real_t sum=0,W;
//  BeamProfile=(DATA_TYPE*)get_spc(BEAM_RESOLUTION,sizeof(DATA_TYPE));
  size_t dims[1] = {m_AdvParams->BEAM_RESOLUTION };
  BeamProfile = RealArrayType::New(dims, "BeamProfile");
  W=BEAM_WIDTH/2;
  for (i=0; i <m_AdvParams->BEAM_RESOLUTION ;i++)
  {
    //BeamProfile->d[i] = (1.0/(BEAM_WIDTH)) * ( 1 + cos ((PI/W)*fabs(-W + i*(BEAM_WIDTH/BEAM_RESOLUTION))));
    BeamProfile->d[i] = 0.54 - 0.46*cos((2.0*M_PI/m_AdvParams->BEAM_RESOLUTION)*i);
    sum=sum+BeamProfile->d[i];
  }

  //Normalize the beam to have an area of 1

  for (i=0; i <m_AdvParams->BEAM_RESOLUTION ;i++)
  {

    BeamProfile->d[i]/=sum;
    BeamProfile->d[i]/=m_Sinogram->delta_t;//This is for proper normalization
    // printf("%lf\n",BeamProfile->d[i]);
  }



}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t SOCEngine::computeCost(RealVolumeType::Pointer ErrorSino,RealVolumeType::Pointer Weight)
{
  Real_t cost=0,temp=0;
  Real_t delta;
  Real_t errSinoValue = 0.0;
#ifdef EIMTOMO_USE_QGGMRF
  //DATA_TYPE MRF_C_TIMES_SIGMA_P_Q= MRF_C*SIGMA_X_P_Q;
#endif
//  int16_t p,q,r;

//Data Mismatch Error
  for (int16_t i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (int16_t j = 0; j < m_Sinogram->N_r; j++)
    {
      for (int16_t k = 0; k < m_Sinogram->N_t; k++)
      {
        errSinoValue = ErrorSino->getValue(i, j, k);
        cost += (errSinoValue * errSinoValue * Weight->getValue(i, j, k) );
      }
    }
  }

  cost /= 2;

//  std::cout << "\nCompute Cost: Data mismatch term = " << cost;
//  fflush(stdout);

//Prior Model Error
  temp=0;
#ifndef EIMTOMO_USE_QGGMRF
  for (int16_t i = 0; i < m_Geometry->N_z; i++)
    for (int16_t j = 0; j < m_Geometry->N_x; j++)
      for(int16_t k = 0; k < m_Geometry->N_y; k++)
      {

        if(k+1 <  m_Geometry->N_y)
          temp += FILTER[2][1][1]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j][k+1]),MRF_P);


        if(j+1 < m_Geometry->N_x)
        {
          if(k-1 >= 0)
            temp += FILTER[0][1][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k-1]),MRF_P);


          temp += FILTER[1][1][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k]),MRF_P);


          if(k+1 < m_Geometry->N_y)
            temp += FILTER[2][1][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k+1]),MRF_P);

        }

        if(i+1 < m_Geometry->N_z)
        {

          if(j-1 >= 0)
            temp += FILTER[1][2][0]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k]),MRF_P);

          temp += FILTER[1][2][1]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k]),MRF_P);

          if(j+1 < m_Geometry->N_x)
            temp += FILTER[1][2][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k]),MRF_P);


          if(j-1 >= 0)
          {
            if(k-1 >= 0)
              temp += FILTER[0][2][0]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k-1]),MRF_P);

            if(k+1 < m_Geometry->N_y)
              temp += FILTER[2][2][0]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k+1]),MRF_P);

          }

          if(k-1 >= 0)
            temp += FILTER[0][2][1]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k-1]),MRF_P);

          if(j+1 < m_Geometry->N_x)
          {
            if(k-1 >= 0)
              temp += FILTER[0][2][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k-1]),MRF_P);

            if(k+1 < m_Geometry->N_y)
              temp+= FILTER[2][2][2]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k+1]),MRF_P);
          }

          if(k+1 < m_Geometry->N_y)
            temp+= FILTER[2][2][1]*pow(fabs(m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k+1]),MRF_P);
        }
      }
    cost+=(temp/(MRF_P*SIGMA_X_P));
#else

  for (int16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (int16_t j = 0; j < m_Geometry->N_x; j++)
    {
      for (int16_t k = 0; k < m_Geometry->N_y; k++)
      {

        if(k + 1 < m_Geometry->N_y)
        {
          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j, k + 1);
          temp += FILTER[INDEX_3(2,1,1)] * QGGMRF::Value(delta, &m_QGGMRF_Values);

        }

        if(j + 1 < m_Geometry->N_x)
        {
          if(k - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k - 1);
            temp += FILTER[INDEX_3(0,1,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }

          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k);
          temp += FILTER[INDEX_3(1,1,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);

          if(k + 1 < m_Geometry->N_y)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k + 1);
            temp += FILTER[INDEX_3(2,1,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }

        }

        if(i + 1 < m_Geometry->N_z)
        {

          if(j - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k);
            temp += FILTER[INDEX_3(1,2,0)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }

          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k);
          temp += FILTER[INDEX_3(1,2,1)] * QGGMRF::Value(delta, &m_QGGMRF_Values);

          if(j + 1 < m_Geometry->N_x)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k);
            temp += FILTER[INDEX_3(1,2,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }

          if(j - 1 >= 0)
          {
            if(k - 1 >= 0)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k - 1);
              temp += FILTER[INDEX_3(0,2,0)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
            }

            if(k + 1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k + 1);
              temp += FILTER[INDEX_3(2,2,0)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
            }

          }

          if(k - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k - 1);
            temp += FILTER[INDEX_3(0,2,1)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }

          if(j + 1 < m_Geometry->N_x)
          {
            if(k - 1 >= 0)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k - 1);
              temp += FILTER[INDEX_3(0,2,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
            }

            if(k + 1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k + 1);
              temp += FILTER[INDEX_3(2,2,2)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
            }
          }

          if(k + 1 < m_Geometry->N_y)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k + 1);
            temp += FILTER[INDEX_3(2,2,1)] * QGGMRF::Value(delta, &m_QGGMRF_Values);
          }
        }
      }
    }
  }
  cost+=(temp);
#endif //QGGMRF

  //printf("Cost calculation End..\n");


//Noise Error
  if(m_AdvParams->NOISE_MODEL) {
  temp = 0;
  for (int16_t i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (int16_t j = 0; j < m_Sinogram->N_r; j++)
    {
      for (int16_t k = 0; k < m_Sinogram->N_t; k++)
      {
        if(Weight->getValue(i, j, k) != 0)
        temp+= log(2 * M_PI * (1.0/Weight->getValue(i, j, k)));
      }
    }
  }
  temp/=2;
  cost += temp;
}//NOISE_MODEL
  return cost;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AMatrixCol::Pointer SOCEngine::calculateAMatrixColumnPartial(uint16_t row,uint16_t col, uint16_t slice,
                                                                 RealVolumeType::Pointer DetectorResponse)
{
  int32_t j, k, sliceidx;
  Real_t x, z, y;
  Real_t r; //this is used to find where does the ray passing through the voxel at certain angle hit the detector
  Real_t t; //this is similar to r but along the y direction
  Real_t tmin, tmax;
  Real_t rmax, rmin; //stores the start and end points of the pixel profile on the detector
  Real_t R_Center, TempConst, checksum = 0, delta_r;
//  DATA_TYPE Integral = 0;
  Real_t T_Center, delta_t;
  Real_t MaximumSpacePerColumn; //we will use this to allocate space
  Real_t AvgNumXElements, AvgNumYElements; //This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
//  DATA_TYPE ProfileThickness,stepsize;

  //interpolation variables
  Real_t w1, w2, w3, w4, f1, InterpolatedValue, ContributionAlongT;
//  DATA_TYPE f2;
  int32_t index_min, index_max, slice_index_min, slice_index_max, index_delta_r, index_delta_t; //stores the detector index in which the profile lies
  int32_t BaseIndex, FinalIndex;
//  int32_t ProfileIndex=0;
//  int32_t NumOfDisplacements=32;
  uint32_t count = 0;

  sliceidx = 0;



  x = m_Geometry->x0 + ((Real_t)col + 0.5) * m_TomoInputs->delta_xz; //0.5 is for center of voxel. x_0 is the left corner
  z = m_Geometry->z0 + ((Real_t)row + 0.5) * m_TomoInputs->delta_xz; //0.5 is for center of voxel. x_0 is the left corner
  y = m_Geometry->y0 + ((Real_t)slice + 0.5) * m_TomoInputs->delta_xy;

  TempConst = (m_AdvParams->PROFILE_RESOLUTION) / (2 * m_TomoInputs->delta_xz);

  //alternately over estimate the maximum size require for a single AMatrix column
  AvgNumXElements = ceil(3 * m_TomoInputs->delta_xz / m_Sinogram->delta_r);
  AvgNumYElements = ceil(3 * m_TomoInputs->delta_xy / m_Sinogram->delta_t);
  MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements) * m_Sinogram->N_theta;

  size_t dims[1] = { MaximumSpacePerColumn };
  AMatrixCol::Pointer Temp = AMatrixCol::New(dims, 0);
//  AMatrixCol* Temp = (AMatrixCol*)get_spc(1, sizeof(AMatrixCol)); //This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end
//
//  Temp->values = (Real_t*)get_spc((uint32_t)MaximumSpacePerColumn, sizeof(Real_t));
//  Temp->index = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn, sizeof(uint32_t));

  if(m_AdvParams->AREA_WEIGHTED)
  {
    for (uint32_t i = 0; i < m_Sinogram->N_theta; i++)
    {

      r = x * cosine->d[i] - z * sine->d[i];
      t = y;

      rmin = r - m_TomoInputs->delta_xz;
      rmax = r + m_TomoInputs->delta_xz;

      tmin = (t - m_TomoInputs->delta_xy / 2) > m_Sinogram->T0 ? t - m_TomoInputs->delta_xy / 2 : m_Sinogram->T0;
      tmax = (t + m_TomoInputs->delta_xy / 2) <= m_Sinogram->TMax ? t + m_TomoInputs->delta_xy / 2 : m_Sinogram->TMax;

      if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax) continue;

      index_min = static_cast<int32_t>(floor(((rmin - m_Sinogram->R0) / m_Sinogram->delta_r)));
      index_max = static_cast<int32_t>(floor((rmax - m_Sinogram->R0) / m_Sinogram->delta_r));

      if(index_max >= m_Sinogram->N_r) index_max = m_Sinogram->N_r - 1;

      if(index_min < 0) index_min = 0;

      slice_index_min = static_cast<int32_t>(floor((tmin - m_Sinogram->T0) / m_Sinogram->delta_t));
      slice_index_max = static_cast<int32_t>(floor((tmax - m_Sinogram->T0) / m_Sinogram->delta_t));

      if(slice_index_min < 0) slice_index_min = 0;
      if(slice_index_max >= m_Sinogram->N_t) slice_index_max = m_Sinogram->N_t - 1;

      BaseIndex = i * m_Sinogram->N_r; //*Sinogram->N_t;

      for (j = index_min; j <= index_max; j++) //Check
      {

        //Accounting for Beam width
        R_Center = (m_Sinogram->R0 + (((Real_t)j) + 0.5) * (m_Sinogram->delta_r)); //the 0.5 is to get to the center of the detector

        //Find the difference between the center of detector and center of projection and compute the Index to look up into
        delta_r = fabs(r - R_Center);
        index_delta_r = static_cast<int32_t>(floor((delta_r / OffsetR)));

        if(index_delta_r >= 0 && index_delta_r < m_AdvParams->DETECTOR_RESPONSE_BINS)
        {
          T_Center = (m_Sinogram->T0 + (((Real_t)sliceidx) + 0.5) * (m_Sinogram->delta_t));
          delta_t = fabs(t - T_Center);
          index_delta_t = 0; //floor(delta_t/OffsetT);
          if(index_delta_t >= 0 && index_delta_t < m_AdvParams->DETECTOR_RESPONSE_BINS)
          {
            //Using index_delta_t,index_delta_t+1,index_delta_r and index_delta_r+1 do bilinear interpolation
            w1 = delta_r - index_delta_r * OffsetR;
            w2 = (index_delta_r + 1) * OffsetR - delta_r;

            w3 = delta_t - index_delta_t * OffsetT;
            w4 = (index_delta_r + 1) * OffsetT - delta_t;

            uint16_t iidx = index_delta_r + 1 < m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_r + 1 : m_AdvParams->DETECTOR_RESPONSE_BINS - 1;
            f1 = (w2 / OffsetR) * DetectorResponse->getValue(index_delta_t, i, index_delta_r)
                + (w1 / OffsetR) * DetectorResponse->getValue(index_delta_t, i, iidx);
            //  f2 = (w2/OffsetR)*DetectorResponse[index_delta_t+1 < m_AdvParams->DETECTOR_RESPONSE_BINS ?index_delta_t+1 : m_AdvParams->DETECTOR_RESPONSE_BINS-1][i][index_delta_r] + (w1/OffsetR)*DetectorResponse[index_delta_t+1 < m_AdvParams->DETECTOR_RESPONSE_BINS? index_delta_t+1:m_AdvParams->DETECTOR_RESPONSE_BINS][i][index_delta_r+1 < m_AdvParams->DETECTOR_RESPONSE_BINS? index_delta_r+1:m_AdvParams->DETECTOR_RESPONSE_BINS-1];

            if(sliceidx == slice_index_min) ContributionAlongT = (sliceidx + 1) * m_Sinogram->delta_t - tmin;
            else if(sliceidx == slice_index_max) ContributionAlongT = tmax - (sliceidx) * m_Sinogram->delta_t;
            else
            {
              ContributionAlongT = m_Sinogram->delta_t;
            }
            InterpolatedValue = f1; //*ContributionAlongT;//(w3/OffsetT)*f2 + (w4/OffsetT)*f2;
            if(InterpolatedValue > 0)
            {
              FinalIndex = BaseIndex + (int32_t)j; //+ (int32_t)sliceidx * Sinogram->N_r;
              Temp->values[count] = InterpolatedValue; //DetectorResponse[index_delta_t][i][index_delta_r];
              Temp->index[count] = FinalIndex; //can instead store a triple (row,col,slice) for the sinogram
              count++;
            }
          }
        }
      }
    }
  }

  //AMatrixCol* Ai = (AMatrixCol*)get_spc(1, sizeof(AMatrixCol));

  dims[0] = count;
  AMatrixCol::Pointer Ai = AMatrixCol::New(dims, 0);
//
//  Ai->values = (Real_t*)get_spc(count, sizeof(Real_t));
//  Ai->index = (uint32_t*)get_spc(count, sizeof(uint32_t));
  k = 0;
  for (uint32_t i = 0; i < count; i++)
  {
    if(Temp->values[i] > 0.0)
    {
      Ai->values[k] = Temp->values[i];
      checksum += Ai->values[k];
      Ai->index[k] = Temp->index[i];
      k++;
    }
  }
  Ai->setCount(k);

//  free(Temp->values);
//  free(Temp->index);
//  free(Temp);
  return Ai;
}

#ifndef EIMTOMO_USE_QGGMRF
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double SOCEngine::surrogateFunctionBasedMin(Real_t currentVoxelValue)
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
Real_t SOCEngine::absMaxArray(std::vector<Real_t> &Array)
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
void SOCEngine::ComputeVSC()
{
  Real_t filter_op = 0;
  int err = 0;
  FILE *Fp = NULL;
  MAKE_OUTPUT_FILE(Fp, err, m_TomoInputs->tempDir, ScaleOffsetCorrection::MagnitudeMapFile);
  if(err < 0)
  {

  }
  fwrite( MagUpdateMap->getPointer(0, 0), m_Geometry->N_x * m_Geometry->N_z, sizeof(Real_t), Fp);
  fclose(Fp);

  // std::cout<<"Starting to filter the magnitude"<<std::endl;
  // std::cout<<m_Geometry->N_x<<" " <<m_Geometry->N_z<<std::endl;
  for (int16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (int16_t j = 0; j < m_Geometry->N_x; j++)
    {
      filter_op = 0;
      for (int16_t p = -2; p <= 2; p++)
      {
        for (int16_t q = -2; q <= 2; q++)
        {
          if(i + p >= 0 && i + p < m_Geometry->N_z && j + q >= 0 && j + q < m_Geometry->N_x)
          {
            filter_op += HAMMING_WINDOW[p + 2][q + 2] * MagUpdateMap->getValue(i + p, j + q);
          }
        }
      }
      FiltMagUpdateMap->setValue(filter_op, i, j);
    }
  }

  for (int16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (int16_t j = 0; j < m_Geometry->N_x; j++)
    {
     //MagUpdateMap->d[i][j]=FiltMagUpdateMap->d[i][j];
      MagUpdateMap->setValue(FiltMagUpdateMap->getValue(i, j), i, j);
    }
  }

  MAKE_OUTPUT_FILE(Fp, err, m_TomoInputs->tempDir, ScaleOffsetCorrection::FilteredMagMapFile);
  if(err < 0)
  {

  }
  fwrite( FiltMagUpdateMap->getPointer(0, 0), m_Geometry->N_x * m_Geometry->N_z, sizeof(Real_t), Fp);
  fclose(Fp);
}


//Sort the entries of FiltMagUpdateMap and set the threshold to be ? percentile
Real_t SOCEngine::SetNonHomThreshold()
{
  size_t dims[2] =
  { m_Geometry->N_z * m_Geometry->N_x, 0 };
  RealArrayType::Pointer TempMagMap = RealArrayType::New(dims, "TempMagMap");

  uint32_t ArrLength = m_Geometry->N_z * m_Geometry->N_x;
  Real_t threshold;

  //Copy into a linear list for easier partial sorting
  for (uint32_t i = 0; i < m_Geometry->N_z; i++)
    for (uint32_t j = 0; j < m_Geometry->N_x; j++)
    {
      //TempMagMap->d[i*m_Geometry->N_x+j]=i*m_Geometry->N_x+j;
      TempMagMap->d[i * (uint32_t)m_Geometry->N_x + j] = MagUpdateMap->getValue(i, j);
    }

  uint16_t percentile_index = ArrLength / NUM_NON_HOMOGENOUS_ITER;
  //Partial selection sort

  Real_t max;
  uint32_t max_index;
  for (uint32_t i = 0; i <= percentile_index; i++)
  {
    max = TempMagMap->d[i];
    max_index = i;
    for (uint32_t j = i + 1; j < ArrLength; j++)
    {
      if(TempMagMap->d[j] > max)
      {
        max = TempMagMap->d[j];
        max_index = j;
      }
    }
    Real_t temp = TempMagMap->d[i];
    TempMagMap->d[i] = TempMagMap->d[max_index];
    TempMagMap->d[max_index] = temp;
  }

  //Insertion sort
  /*
   int32_t j;
   DATA_TYPE key;
   for (uint32_t i=1 ; i < ArrLength; i++)
   {
   j=i-1;
   key=TempMagMap->d[i];
   while( j>=0 && TempMagMap->d[j] < key)
   {
   TempMagMap->d[j+1]=TempMagMap->d[j];
   j--;
   }
   TempMagMap->d[j+1]=key;
   }

   //TempMagMap is a local variable and will clean up its own memory when this method exits
   uint16_t percentile_index=ArrLength/NUM_NON_HOMOGENOUS_ITER;
   std::cout<<ArrLength<<" "<<percentile_index<<std::endl;
   */
  threshold = TempMagMap->d[percentile_index];
  return threshold;
}




