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

#if TomoEngine_USE_PARALLEL_ALGORITHMS
//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>
//#include <tbb/atomic.h>
//#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif

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
#include "TomoEngine/mt/mt19937ar.h"
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
   v->excludedViews;
   v->goodViews;
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
void SOCEngine::InitializeScaleOffsetParams(ScaleOffsetParamsPtr v)
{
  v->I_0 = RealArrayType::NullPointer();
  v->mu = RealArrayType::NullPointer();
  v->alpha = RealArrayType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SOCEngine::SOCEngine()
{
  initVariables();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SOCEngine::~SOCEngine()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::initVariables()
{
  FILTER[0][0][0] = 0.0302; FILTER[0][0][1] = 0.0370; FILTER[0][0][2] = 0.0302;
  FILTER[0][1][0] = 0.0370; FILTER[0][1][1] = 0.0523; FILTER[0][1][2] = 0.0370;
  FILTER[0][2][0] = 0.0302; FILTER[0][2][1] = 0.0370; FILTER[0][2][2] = 0.0302;

  FILTER[1][0][0] = 0.0370; FILTER[1][0][1] = 0.0523; FILTER[1][0][2] = 0.0370;
  FILTER[1][1][0] = 0.0523; FILTER[1][1][1] = 0.0000; FILTER[1][1][2] = 0.0523;
  FILTER[1][2][0] = 0.0370; FILTER[1][2][1] = 0.0523; FILTER[1][2][2] = 0.0370;

  FILTER[2][0][0] = 0.0302; FILTER[2][0][1] = 0.0370; FILTER[2][0][2] = 0.0302;
  FILTER[2][1][0] = 0.0370; FILTER[2][1][1] = 0.0523; FILTER[2][1][2] = 0.0370;
  FILTER[2][2][0] = 0.0302; FILTER[2][2][1] = 0.0370; FILTER[2][2][2] = 0.0302;

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

 // int16_t i,j,k,Idx;
  size_t dims[3];

  Int32ArrayType::Pointer Counter;
  UInt8Image_t::Pointer VisitCount;


  uint16_t MaxNumberOfDetectorElts;
  uint8_t status;//set to 1 if ICD has converged


  Real_t checksum = 0,temp;

  RealImage_t::Pointer VoxelProfile;
  RealVolumeType::Pointer detectorResponse;
  RealVolumeType::Pointer H_t;


  RealVolumeType::Pointer Y_Est;//Estimated Sinogram
  RealVolumeType::Pointer Final_Sinogram;//To store and write the final sinogram resulting from our reconstruction
  RealVolumeType::Pointer ErrorSino;//Error Sinogram
  RealVolumeType::Pointer Weight;//This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation

  RNGVars* RandomNumber;
  std::string indent("");

#ifdef COST_CALCULATE
  std::string filepath(m_TomoInputs->tempDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::CostFunctionFile);

  CostData::Pointer cost = CostData::New();
  cost->initOutputFile(filepath);
#endif


#if TomoEngine_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
#endif

  // Initialize the Sinogram
  if (m_TomoInputs == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The TomoInput Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }
  //Based on the inputs , calculate the "other" variables in the structure definition
  if (m_Sinogram == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The Sinogram Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }
  if (m_Geometry == NULL)
  {
    setErrorCondition(-1);
    notify("Error: The Geometry Structure was NULL. The proper API is to supply this class with that structure,", 100, Observable::UpdateProgressValueAndMessage);
    return;
  }


  // Read the Input data from the supplied data file
  err = readInputData();
  if (err < 0)
  {
    return;
  }

  err = initializeBrightFieldData();
  if (err < 0)
  {
    return;
  }

  err = createInitialGainsData();
  if (err < 0)
  {
    return;
  }
  err = createInitialOffsetsData();
  if (err < 0)
  {
    return;
  }
  err = createInitialVariancesData();
  if (err < 0)
  {
    return;
  }

#ifdef DEBUG
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
#endif

  // Initialize the Geometry data from a rough reconstruction
  err = initializeRoughReconstructionData();
  if (err < 0)
  {
    return;
  }

  //Gain, Offset and Variance Parameter Structures
  ScaleOffsetParamsPtr NuisanceParams = ScaleOffsetParamsPtr(new ScaleOffsetParams);
  dims[1] = m_Sinogram->N_t;
  dims[0] = m_Sinogram->N_theta;
  NuisanceParams->I_0 = RealArrayType::New(dims, "NuisanceParams->I_0");
  NuisanceParams->mu = RealArrayType::New(dims, "NuisanceParams->mu");
#ifdef NOISE_MODEL
  //alpha is the noise variance adjustment factor
  NuisanceParams->alpha = RealArrayType::New(dims, "NuisanceParams->alpha");
#else
  NuisanceParams->alpha = RealArrayType::NullPointer();
#endif

  // initialize variables
  uint16_t Idx = 0;

#ifdef WRITE_INTERMEDIATE_RESULTS
  Real_t Fraction = 0.1;//write this fraction of the iterations
  int16_t NumOfWrites = floor((Real_t)(m_TomoInputs->NumIter)*Fraction);
  int16_t WriteCount = 0;
  char Filename[100];
  char buffer[20];
  void* TempPointer;
  size_t NumOfBytesWritten;
#endif

#ifdef ROI
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
  OffsetR = ((m_TomoInputs->delta_xz/sqrt(3.0)) + m_Sinogram->delta_r/2)/DETECTOR_RESPONSE_BINS;
  OffsetT = ((m_TomoInputs->delta_xz/2) + m_Sinogram->delta_t/2)/DETECTOR_RESPONSE_BINS;

  BEAM_WIDTH = m_Sinogram->delta_r;


#ifndef QGGMRF
  MRF_P = m_TomoInputs->p;
  SIGMA_X_P = pow(m_TomoInputs->SigmaX,MRF_P);
#else
  MRF_P = 2;
  MRF_Q = m_TomoInputs->p;
  MRF_C = 0.01;
  MRF_ALPHA = 1.5;
  SIGMA_X_P = pow(m_TomoInputs->SigmaX,MRF_P);
  SIGMA_X_P_Q = pow(m_TomoInputs->SigmaX, (MRF_P - MRF_Q));
  SIGMA_X_Q = pow(m_TomoInputs->SigmaX,MRF_Q);
#endif //QGGMRF

  //globals assosiated with finding the optimal gain and offset parameters

  dims[0] = m_Sinogram->N_theta;
  dims[1] = 3;
  dims[2] = 0;
  //Hold the coefficients of a quadratic equation
  QuadraticParameters = RealImage_t::New(dims, "QuadraticParameters");
  Qk_cost = RealImage_t::New(dims, "Qk_cost");
  dims[1] = 2;
  bk_cost = RealImage_t::New(dims, "bk_cost");

  dims[0] = m_Sinogram->N_theta;
  ck_cost = RealArrayType::New(dims, "ck_cost");

  NumOfViews = m_Sinogram->N_theta;
  LogGain = m_Sinogram->N_theta*log(m_Sinogram->targetGain);
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
  dResponseFilter->setBeamWidth(BEAM_WIDTH);
  dResponseFilter->setOffsetR(OffsetR);
  dResponseFilter->setOffsetT(OffsetT);
  dResponseFilter->setVoxelProfile(VoxelProfile);
  dResponseFilter->setBeamProfile(BeamProfile);
  dResponseFilter->setObservers(getObservers());
  dResponseFilter->execute();
  if (dResponseFilter->getErrorCondition() < 0)
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
  responseWriter->setObservers(getObservers());
  responseWriter->setResponse(detectorResponse);
  responseWriter->execute();
  if (responseWriter->getErrorCondition() < 0)
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
  ::memset( VisitCount->d, 0, dims[0] * dims[1] * sizeof(uint8_t));
#endif//Random update


  dims[0]=m_Geometry->N_z;//height
  dims[1]=m_Geometry->N_x;//width
  dims[2]=0;

  MagUpdateMap = RealImage_t::New(dims, "Update Map for voxel lines");
  FiltMagUpdateMap = RealImage_t::New(dims, "Update Map for voxel lines");
  MagUpdateMask = UInt8Image_t::New(dims, "Update Mask for selecting voxel lines NHICD");

#ifdef ROI
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
  dims[2] = DETECTOR_RESPONSE_BINS;
  H_t = RealVolumeType::New(dims, "H_t");
  initializeHt(H_t);


  checksum=0;

#ifdef STORE_A_MATRIX

  AMatrixCol**** AMatrix = (AMatrixCol ****)multialloc(sizeof(AMatrixCol*),3,m_Geometry->N_y,m_Geometry->N_z,m_Geometry->N_x);
#else

  //TODO: All this needs to be deallocated at some point
//  DATA_TYPE y;
//  DATA_TYPE t, tmin, tmax, ProfileThickness;
//  int16_t slice_index_min, slice_index_max;
  AMatrixCol** TempCol = (AMatrixCol**)get_spc(m_Geometry->N_x * m_Geometry->N_z, sizeof(AMatrixCol*));
  AMatrixCol* VoxelLineResponse = (AMatrixCol*)get_spc(m_Geometry->N_y, sizeof(AMatrixCol));
  MaxNumberOfDetectorElts = (uint16_t)((m_TomoInputs->delta_xy / m_Sinogram->delta_t) + 2);
  for (uint16_t i = 0; i < m_Geometry->N_y; i++)
  {
    VoxelLineResponse[i].count = 0;
    VoxelLineResponse[i].values = (Real_t*)get_spc(MaxNumberOfDetectorElts, sizeof(Real_t));
    VoxelLineResponse[i].index = (uint32_t*)get_spc(MaxNumberOfDetectorElts, sizeof(uint32_t));
  }
#endif

  //Calculating A-Matrix one column at a time
  //For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
  // And then store only the non zero entries by allocating a new array of the desired size
  //k=0;

  checksum = 0;
  //q = 0;

#ifdef STORE_A_MATRIX
  for(uint16_t i = 0;i < m_Geometry->N_y; i++)
  {
    for(uint16_t j = 0;j < m_Geometry->N_z; j++)
    {
      for (uint16_t k = 0; k < m_Geometry->N_x; k++)
      {
        //  AMatrix[q++]=CE_CalculateAMatrixColumn(i,j,Sinogram,Geometry,VoxelProfile);
        AMatrix[i][j][k]=calculateAMatrixColumn(j,k,i,m_Sinogram,m_Geometry,VoxelProfile);//row,col,slice

        for(p = 0; p < AMatrix[i][j][k]->count; p++)
        {
          checksum += AMatrix[i][j][k]->values[p];}
        //   printf("(%d,%d,%d) %lf \n",i,j,k,AMatrix[i][j][k]->values);
        checksum = 0;
      }}}
  printf("Stored A matrix\n");
#else
  temp = 0;
  uint32_t voxel_count = 0;
  for (uint16_t z = 0; z < m_Geometry->N_z; z++)
  {
    for (uint16_t x = 0; x < m_Geometry->N_x; x++)
    {
      TempCol[voxel_count] = (AMatrixCol*)calculateAMatrixColumnPartial(z, x, 0, detectorResponse);
      temp += TempCol[voxel_count]->count;
      if(0 == TempCol[voxel_count]->count)
      {
        //If this line is never hit and the Object is badly initialized
        //set it to zero
        for (uint16_t y = 0; y < m_Geometry->N_y; y++)
        {
          //m_Geometry->Object->d[z][x][y] = 0;
          m_Geometry->Object->setValue(0.0, z, x, y);
        }
      }
      voxel_count++;
    }
  }
#endif


  storeVoxelResponse(H_t, VoxelLineResponse);

  printf("Number of non zero entries of the forward projector is %lf\n",temp);
  printf("Geometry-Z %d\n",m_Geometry->N_z);


  //Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
  //is Y_Est initailized to zero?
  initializeVolume(Y_Est, 0.0);

  RandomNumber=init_genrand(1ul);

  notify("Starting Forward Projection", 10, Observable::UpdateProgressValueAndMessage);
  START_TIMER;
  // This next section looks crazy with all the #if's but this makes sure we are
  // running the exact same code whether in parallel or serial.

#if TomoEngine_USE_PARALLEL_ALGORITHMS
  tbb::task_group* g = new tbb::task_group;
  std::cout << "Default Number of Threads to Use: " << init.default_num_threads() << std::endl;
  std::cout << "Forward Projection Running in Parallel." << std::endl;
#else
  std::cout << "Forward Projection Running in Serial." << std::endl;
#endif
  // Queue up a thread for each z layer of the Geometry. The threads will only be
  // run as hardware resources open up so this will not just fire up a gazillion
  // threads.
  for (uint16_t t = 0; t < m_Geometry->N_z; t++)
  {
#if TomoEngine_USE_PARALLEL_ALGORITHMS
    g->run(ForwardProject(m_Sinogram, m_Geometry, TempCol,VoxelLineResponse,Y_Est, &NuisanceParams, t));
#else
    ForwardProject fp(m_Sinogram.get(), m_Geometry.get(), TempCol,VoxelLineResponse,Y_Est, NuisanceParams.get(), t);
    fp();
#endif

  }
#if TomoEngine_USE_PARALLEL_ALGORITHMS
  g->wait(); // Wait for all the threads to complete before moving on.
  delete g;
#endif
  printf("\n");
  STOP_TIMER;
  PRINT_TIME("Forward Project Time");


  //Calculate Error Sinogram - Can this be combined with previous loop?
  //Also compute weights of the diagonal covariance matrix
  calculateMeasurementWeight(Weight, NuisanceParams, ErrorSino, Y_Est);


#ifdef FORWARD_PROJECT_MODE
  return 0;//exit the program once we finish forward projecting the object
#endif//Forward Project mode


#ifdef COST_CALCULATE
  err = calculateCost(cost, Weight, ErrorSino);
#endif //Cost calculation endif

//  int totalLoops = m_TomoInputs->NumOuterIter * m_TomoInputs->NumIter;
  std::stringstream ss;
  //Loop through every voxel updating it by solving a cost function
  for(int16_t OuterIter = 0; OuterIter < m_TomoInputs->NumOuterIter; OuterIter++)
  {
    ss.str(""); // Clear the string stream
    indent = "";

	//The first time we may need to update voxels multiple times and then on just optimize over I,d,\sigma,f once each outer loop
    if(OuterIter != 0)
    {
      m_TomoInputs->NumIter = 1;
    }



    for (int16_t Iter = 0; Iter < m_TomoInputs->NumIter; Iter++)
    {

	  std::cout<<OuterIter<<"/"<<m_TomoInputs->NumOuterIter<<" "<<Iter<<"/"<<m_TomoInputs->NumIter<<std::endl;
      indent = "  ";
//      ss << "Outer Iteration: " << OuterIter << " of " << m_TomoInputs->NumOuterIter;
//      ss << "   Inner Iteration: " << Iter << " of " << m_TomoInputs->NumIter;
//      float currentLoop = OuterIter * m_TomoInputs->NumIter + Iter;
//      notify(ss.str(), currentLoop / totalLoops * 100);
      indent = "    ";
      // This is all done PRIOR to calling what will become a method
      VoxelUpdateType updateType = RegularRandomOrderUpdate;
#ifdef NHICD
      if(0 == Iter % 2)
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
       status =
          updateVoxels(OuterIter, Iter, updateType, VisitCount, RandomNumber,
                       TempCol, ErrorSino, Weight, VoxelLineResponse, NuisanceParams.get(), Mask, cost);

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

    } /* ++++++++++ END Inner Iteration Loop +++++++++++++++ */
#ifdef JOINT_ESTIMATION
     err = jointEstimation(Weight, NuisanceParams, ErrorSino, Y_Est, cost);
     if (err < 0)
     {
       break;
     }
#endif//Joint estimation endif

#ifdef NOISE_MODEL
    updateWeights(Weight, NuisanceParams, ErrorSino);
#ifdef COST_CALCULATE
    err = calculateCost(cost, Weight, ErrorSino);
    if (err < 0)
    {
	  std::cout<<"Cost went up after variance update"<<std::endl;
      break;
    }
#endif//cost
    if(0 == status && OuterIter >= 1)//&& VarRatio < STOPPING_THRESHOLD_Var_k && I_kRatio < STOPPING_THRESHOLD_I_k && Delta_kRatio < STOPPING_THRESHOLD_Delta_k)
	{
	 std::cout<<"Exiting the code because status =0"<<std::endl;
      break;
	}
#else
    if(0 == status)//&& I_kRatio < STOPPING_THRESHOLD_I_k && Delta_kRatio < STOPPING_THRESHOLD_Delta_k)
      break;
#endif //Noise Model

  }/* ++++++++++ END Outer Iteration Loop +++++++++++++++ */

  indent = "";
#ifdef DEBUG
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


#ifdef DEBUG_CONSTRAINT_OPT
  FILE *Fp8 = NULL;
  int fileError=0;
  MAKE_OUTPUT_FILE(Fp8, fileError, m_TomoInputs->tempDir, ScaleOffsetCorrection::CostFunctionCoefficientsFile);
  if (fileError >= 0)
  {
    // Write the output File
  }
#endif


/* Write the Gains and Offsets to an output file */
  writeNuisanceParameters(NuisanceParams);

#ifdef DEBUG
  std::cout << "Tilt\tFinal Gains\tFinal Offsets\tFinal Variances" << std::endl;
  for (uint16_t i_theta = 0; i_theta < getSinogram()->N_theta; i_theta++)
  {
    std::cout << i_theta << "\t" << NuisanceParams->I_0->d[i_theta] <<
        "\t" << NuisanceParams->mu->d[i_theta] <<
        "\t" << NuisanceParams->alpha->d[i_theta] << std::endl;
  }
#endif


 //calculates Ax and returns a pointer to the memory block
 // Final_Sinogram=forwardProject(detectorResponse, H_t);
  Real_t temp_final = 0.0;
  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
  {
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
//        temp_final = m_Sinogram->counts->d[i_theta][i_r][i_t] - ErrorSino->d[i_theta][i_r][i_t];
        temp_final = m_Sinogram->counts->getValue(i_theta, i_r, i_t) - ErrorSino->getValue(i_theta, i_r, i_t);
        //Final_Sinogram->d[i_theta][i_r][i_t] = temp;
        Final_Sinogram->setValue(temp_final, i_theta, i_r, i_t);
      }
    }
  }



  writeSinogramFile(NuisanceParams, Final_Sinogram);
  writeReconstructionFile();
  writeVtkFile();
  writeMRCFile();



  std::cout << "Final Dimensions of Object: " << std::endl;
  std::cout << "  Nx = " << m_Geometry->N_x << std::endl;
  std::cout << "  Ny = " << m_Geometry->N_y << std::endl;
  std::cout << "  Nz = " << m_Geometry->N_z << std::endl;

  //free(AMatrix);
#ifdef STORE_A_MATRIX
  multifree(AMatrix,2);
  //#else
  //  free((void*)TempCol);
#endif


  notify("Reconstruction Complete", 100, Observable::UpdateProgressValueAndMessage);
  setErrorCondition(0);
  std::cout << "Total Running Time for Execute: " << (EIMTOMO_getMilliSeconds()-totalTime)/1000 << std::endl;
  return;
}


/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void SOCEngine::minMax(Real_t *low,Real_t *high)
{
  uint8_t i,j,k;
  *low=NEIGHBORHOOD[0][0][0];
  *high=NEIGHBORHOOD[0][0][0];

  for(i = 0; i < 3;i++)
  {
    for(j=0; j < 3; j++)
    {
      for(k = 0; k < 3; k++)
      {
        //  if(NEIGHBORHOOD[i][j][k] != 0)
        //  printf("%lf ", NEIGHBORHOOD[i][j][k]);

        if(NEIGHBORHOOD[i][j][k] < *low)
          *low = NEIGHBORHOOD[i][j][k];
        if(NEIGHBORHOOD[i][j][k] > *high)
          *high=NEIGHBORHOOD[i][j][k];
      }
      //  printf("\n");
    }
  }


  if(THETA2 !=0)
  {
  *low = (*low > (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *low);

  *high = (*high < (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *high);
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RealImage_t::Pointer SOCEngine::calculateVoxelProfile()
{
  Real_t angle,MaxValLineIntegral;
  Real_t temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;
  size_t dims[2] = {m_Sinogram->N_theta , PROFILE_RESOLUTION};
  RealImage_t::Pointer VoxProfile = RealImage_t::New(dims, "VoxelProfile");

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

    for(j = 0;j<PROFILE_RESOLUTION;j++)
    {
      t = 2.0*j / PROFILE_RESOLUTION;//2 is the normalized length of the profile (basically equl to 2*delta_xz)
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

            if(index_delta_r < DETECTOR_RESPONSE_BINS)
            {
              w1 = delta_r - index_delta_r * OffsetR;
              w2 = (index_delta_r + 1) * OffsetR - delta_r;
              uint16_t iidx = index_delta_r + 1 < DETECTOR_RESPONSE_BINS ? index_delta_r + 1 : DETECTOR_RESPONSE_BINS - 1;
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

              if(index_delta_t < DETECTOR_RESPONSE_BINS)
              {
                w1 = delta_t - index_delta_t * OffsetT;
                w2 = (index_delta_t + 1) * OffsetT - delta_t;
                uint16_t iidx = index_delta_t + 1 < DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : DETECTOR_RESPONSE_BINS - 1;
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
  size_t dims[1] = { BEAM_RESOLUTION };
  BeamProfile = RealArrayType::New(dims, "BeamProfile");
  W=BEAM_WIDTH/2;
  for (i=0; i < BEAM_RESOLUTION ;i++)
  {
    //BeamProfile->d[i] = (1.0/(BEAM_WIDTH)) * ( 1 + cos ((PI/W)*fabs(-W + i*(BEAM_WIDTH/BEAM_RESOLUTION))));
    BeamProfile->d[i] = 0.54 - 0.46*cos((2*M_PI/BEAM_RESOLUTION)*i);
    sum=sum+BeamProfile->d[i];
  }

  //Normalize the beam to have an area of 1

  for (i=0; i < BEAM_RESOLUTION ;i++)
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
#ifdef QGGMRF
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
#ifndef QGGMRF
  for (i = 0; i < m_Geometry->N_z; i++)
    for (j = 0; j < m_Geometry->N_x; j++)
      for(k = 0; k < m_Geometry->N_y; k++)
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
  /*for (i = 0; i < m_Geometry->N_z; i++)
    for (j = 0; j < m_Geometry->N_x; j++)
      for(k = 0; k < m_Geometry->N_y; k++)
      {

        if(k+1 <  m_Geometry->N_y)
        {
          delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j][k+1];
          temp += FILTER[2][1][1]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q+pow(fabs(delta),MRF_P-MRF_Q));

        }



        if(j+1 < m_Geometry->N_x)
        {
          if(k-1 >= 0)
          {
            delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k-1];
            temp += FILTER[0][1][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q+pow(fabs(delta), MRF_P-MRF_Q));
          }

          delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k];
          temp += FILTER[1][1][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));


          if(k+1 < m_Geometry->N_y)
          {
            delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k+1];
            temp += FILTER[2][1][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
          }

        }

        if(i+1 < m_Geometry->N_z)
        {

          if(j-1 >= 0)
          {
            delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k];
            temp += FILTER[1][2][0]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
          }

          delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k];
          temp += FILTER[1][2][1]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q+pow(fabs(delta), MRF_P-MRF_Q));

          if(j+1 < m_Geometry->N_x)
          {
            delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k];
            temp += FILTER[1][2][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
          }


          if(j-1 >= 0)
          {
            if(k-1 >= 0)
            {
              delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k-1];
              temp += FILTER[0][2][0]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q+pow(fabs(delta), MRF_P-MRF_Q));
            }

            if(k+1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k+1];
              temp += FILTER[2][2][0]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
            }

          }

          if(k-1 >= 0)
          {
            delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k-1];
            temp += FILTER[0][2][1]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
          }

          if(j+1 < m_Geometry->N_x)
          {
            if(k-1 >= 0)
            {
              delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k-1];
              temp += FILTER[0][2][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
            }

            if(k+1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k+1];
              temp+= FILTER[2][2][2]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q + pow(fabs(delta), MRF_P-MRF_Q));
            }
          }

          if(k+1 < m_Geometry->N_y)
          {
            delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k+1];
            temp+= FILTER[2][2][1]*(pow(fabs(delta),MRF_P))/(MRF_C_TIMES_SIGMA_P_Q +pow(fabs(delta), MRF_P-MRF_Q));
          }
        }
      }
      cost+=(temp/SIGMA_X_Q);*/
  for (int16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (int16_t j = 0; j < m_Geometry->N_x; j++)
    {
      for (int16_t k = 0; k < m_Geometry->N_y; k++)
      {

        if(k + 1 < m_Geometry->N_y)
        {
          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j, k + 1);
          temp += FILTER[2][1][1] * CE_QGGMRF_Value(delta);

        }

        if(j + 1 < m_Geometry->N_x)
        {
          if(k - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k - 1);
            temp += FILTER[0][1][2] * CE_QGGMRF_Value(delta);
          }

          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k);
          temp += FILTER[1][1][2] * CE_QGGMRF_Value(delta);

          if(k + 1 < m_Geometry->N_y)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i, j + 1, k + 1);
            temp += FILTER[2][1][2] * CE_QGGMRF_Value(delta);
          }

        }

        if(i + 1 < m_Geometry->N_z)
        {

          if(j - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k);
            temp += FILTER[1][2][0] * CE_QGGMRF_Value(delta);
          }

          delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k);
          temp += FILTER[1][2][1] * CE_QGGMRF_Value(delta);

          if(j + 1 < m_Geometry->N_x)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k);
            temp += FILTER[1][2][2] * CE_QGGMRF_Value(delta);
          }

          if(j - 1 >= 0)
          {
            if(k - 1 >= 0)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k - 1);
              temp += FILTER[0][2][0] * CE_QGGMRF_Value(delta);
            }

            if(k + 1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j - 1, k + 1);
              temp += FILTER[2][2][0] * CE_QGGMRF_Value(delta);
            }

          }

          if(k - 1 >= 0)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k - 1);
            temp += FILTER[0][2][1] * CE_QGGMRF_Value(delta);
          }

          if(j + 1 < m_Geometry->N_x)
          {
            if(k - 1 >= 0)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k - 1);
              temp += FILTER[0][2][2] * CE_QGGMRF_Value(delta);
            }

            if(k + 1 < m_Geometry->N_y)
            {
              delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j + 1, k + 1);
              temp += FILTER[2][2][2] * CE_QGGMRF_Value(delta);
            }
          }

          if(k + 1 < m_Geometry->N_y)
          {
            delta = m_Geometry->Object->getValue(i, j, k) - m_Geometry->Object->getValue(i + 1, j, k + 1);
            temp += FILTER[2][2][1] * CE_QGGMRF_Value(delta);
          }
        }
      }
    }
  }
  cost+=(temp);
#endif //QGGMRF

  //printf("Cost calculation End..\n");


//Noise Error
#ifdef NOISE_MODEL
  temp = 0;
 /* for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    if(Weight->d[i][0][0] != 0)
    {
      temp += log(2 * M_PI * (1.0 / Weight->d[i][0][0])); //2*pi*sigma_k^{2}
    }
  }
  temp *= ((m_Sinogram->N_r * m_Sinogram->N_t) / 2);*/
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
#endif//noise model
  return cost;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void* SOCEngine::calculateAMatrixColumnPartial(uint16_t row,uint16_t col, uint16_t slice,
                                               RealVolumeType::Pointer DetectorResponse)
{
  int32_t j,k,sliceidx;
  Real_t x,z,y;
  Real_t r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
  Real_t t; //this is similar to r but along the y direction
  Real_t tmin,tmax;
  Real_t rmax,rmin;//stores the start and end points of the pixel profile on the detector
  Real_t R_Center,TempConst,checksum = 0,delta_r;
//  DATA_TYPE Integral = 0;
  Real_t T_Center,delta_t;
  Real_t MaximumSpacePerColumn;//we will use this to allocate space
  Real_t AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
//  DATA_TYPE ProfileThickness,stepsize;

  //interpolation variables
  Real_t w1,w2,w3,w4,f1,InterpolatedValue,ContributionAlongT;
//  DATA_TYPE f2;
  int32_t index_min,index_max,slice_index_min,slice_index_max,index_delta_r,index_delta_t;//stores the detector index in which the profile lies
  int32_t BaseIndex,FinalIndex;
//  int32_t ProfileIndex=0;
//  int32_t NumOfDisplacements=32;
  uint32_t count = 0;

  sliceidx = 0;

  AMatrixCol* Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
  AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

  x = m_Geometry->x0 + ((Real_t)col+0.5)*m_TomoInputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  z = m_Geometry->z0 + ((Real_t)row+0.5)*m_TomoInputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
  y = m_Geometry->y0 + ((Real_t)slice + 0.5)*m_TomoInputs->delta_xy;

  TempConst=(PROFILE_RESOLUTION)/(2*m_TomoInputs->delta_xz);

  //alternately over estimate the maximum size require for a single AMatrix column
  AvgNumXElements = ceil(3*m_TomoInputs->delta_xz/m_Sinogram->delta_r);
  AvgNumYElements = ceil(3*m_TomoInputs->delta_xy/m_Sinogram->delta_t);
  MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*m_Sinogram->N_theta;

  Temp->values = (Real_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(Real_t));
  Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));


#ifdef AREA_WEIGHTED
  for(uint32_t i=0;i<m_Sinogram->N_theta;i++)
  {

    r = x*cosine->d[i] - z*sine->d[i];
    t = y;

    rmin = r - m_TomoInputs->delta_xz;
    rmax = r + m_TomoInputs->delta_xz;

    tmin = (t - m_TomoInputs->delta_xy/2) > m_Sinogram->T0 ? t-m_TomoInputs->delta_xy/2 : m_Sinogram->T0;
    tmax = (t + m_TomoInputs->delta_xy/2) <= m_Sinogram->TMax ? t + m_TomoInputs->delta_xy/2 : m_Sinogram->TMax;

    if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax)
      continue;



    index_min = static_cast<int32_t>(floor(((rmin - m_Sinogram->R0)/m_Sinogram->delta_r)));
    index_max = static_cast<int32_t>(floor((rmax - m_Sinogram->R0)/m_Sinogram->delta_r));


    if(index_max >= m_Sinogram->N_r)
      index_max = m_Sinogram->N_r - 1;

    if(index_min < 0)
      index_min = 0;

    slice_index_min = static_cast<int32_t>(floor((tmin - m_Sinogram->T0)/m_Sinogram->delta_t));
    slice_index_max = static_cast<int32_t>(floor((tmax - m_Sinogram->T0)/m_Sinogram->delta_t));

    if(slice_index_min < 0)
      slice_index_min = 0;
    if(slice_index_max >= m_Sinogram->N_t)
      slice_index_max = m_Sinogram->N_t -1;

    BaseIndex = i*m_Sinogram->N_r;//*Sinogram->N_t;

    for(j = index_min;j <= index_max; j++)//Check
    {

      //Accounting for Beam width
      R_Center = (m_Sinogram->R0 + (((Real_t)j) + 0.5) *(m_Sinogram->delta_r));//the 0.5 is to get to the center of the detector

      //Find the difference between the center of detector and center of projection and compute the Index to look up into
      delta_r = fabs(r - R_Center);
      index_delta_r = static_cast<int32_t>(floor((delta_r/OffsetR)));


      if (index_delta_r >= 0 && index_delta_r < DETECTOR_RESPONSE_BINS)
      {

    //    for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
    //    {
          T_Center = (m_Sinogram->T0 + (((Real_t)sliceidx) + 0.5) *(m_Sinogram->delta_t));
          delta_t = fabs(t - T_Center);
          index_delta_t = 0;//floor(delta_t/OffsetT);



          if (index_delta_t >= 0 && index_delta_t < DETECTOR_RESPONSE_BINS)
          {

            //Using index_delta_t,index_delta_t+1,index_delta_r and index_delta_r+1 do bilinear interpolation
            w1 = delta_r - index_delta_r*OffsetR;
            w2 = (index_delta_r+1)*OffsetR - delta_r;

            w3 = delta_t - index_delta_t*OffsetT;
            w4 = (index_delta_r+1)*OffsetT - delta_t;

            uint16_t iidx = index_delta_r+1 < DETECTOR_RESPONSE_BINS ? index_delta_r+1:DETECTOR_RESPONSE_BINS-1;
            f1 = (w2/OffsetR)*DetectorResponse->getValue(index_delta_t, i, index_delta_r)
               + (w1/OffsetR)*DetectorResponse->getValue(index_delta_t, i, iidx);
            //  f2 = (w2/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS ?index_delta_t+1 : DETECTOR_RESPONSE_BINS-1][i][index_delta_r] + (w1/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS? index_delta_t+1:DETECTOR_RESPONSE_BINS][i][index_delta_r+1 < DETECTOR_RESPONSE_BINS? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];

            if(sliceidx == slice_index_min)
              ContributionAlongT = (sliceidx + 1)*m_Sinogram->delta_t - tmin;
            else if(sliceidx == slice_index_max)
              ContributionAlongT = tmax - (sliceidx)*m_Sinogram->delta_t;
            else {
              ContributionAlongT = m_Sinogram->delta_t;
            }


              InterpolatedValue = f1;//*ContributionAlongT;//(w3/OffsetT)*f2 + (w4/OffsetT)*f2;
            if(InterpolatedValue > 0)
            {
              FinalIndex = BaseIndex + (int32_t)j ;//+ (int32_t)sliceidx * Sinogram->N_r;
              Temp->values[count] = InterpolatedValue;//DetectorResponse[index_delta_t][i][index_delta_r];
              Temp->index[count] = FinalIndex;//can instead store a triple (row,col,slice) for the sinogram
              count++;
            }
          }


        //}
      }


    }



  }

#endif


  Ai->values=(Real_t*)get_spc(count,sizeof(Real_t));
  Ai->index=(uint32_t*)get_spc(count,sizeof(uint32_t));
  k=0;
  for(uint32_t i = 0; i < count; i++)
  {
    if(Temp->values[i] > 0.0)
    {
      Ai->values[k]=Temp->values[i];
      checksum+=Ai->values[k];
      Ai->index[k++]=Temp->index[i];
    }

  }

  Ai->count=k;


  free(Temp->values);
  free(Temp->index);
  free(Temp);
  return Ai;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double SOCEngine::surrogateFunctionBasedMin()
{
  double numerator_sum=0;
  double denominator_sum=0;
  double alpha,update=0;
  double product=1;
  uint8_t i,j,k;

  for (i = 0;i < 3;i++)
    for (j = 0; j <3; j++)
      for (k = 0;k < 3; k++)
      {
        if( V != NEIGHBORHOOD[i][j][k])
        {
        product=((double)FILTER[i][j][k]*pow(fabs(V-NEIGHBORHOOD[i][j][k]),(MRF_P-2.0)));
        numerator_sum+=(product*(V-NEIGHBORHOOD[i][j][k]));
        denominator_sum+=product;
        }
      }
  numerator_sum/=SIGMA_X_P;
  denominator_sum/=SIGMA_X_P;

  numerator_sum+=THETA1;
  denominator_sum+=THETA2;

  if(THETA2 > 0)
  {
    alpha=(-1*numerator_sum)/(denominator_sum);
    update = V+Clip(alpha, -V, std::numeric_limits<float>::infinity() );
  }
  else
  {
    update=0;
  }

  if(update > 70000)
    printf("%lf\n",update);

  return update;

}


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
  size_t dims[2]={m_Geometry->N_z*m_Geometry->N_x,0};
  RealArrayType::Pointer TempMagMap=RealArrayType::New(dims, "TempMagMap");

  uint32_t ArrLength=m_Geometry->N_z*m_Geometry->N_x;
  Real_t threshold;

  //Copy into a linear list for easier partial sorting
  for (uint32_t i=0; i < m_Geometry->N_z; i++)
      for (uint32_t j=0; j < m_Geometry->N_x; j++)
    {
      //TempMagMap->d[i*m_Geometry->N_x+j]=i*m_Geometry->N_x+j;
      TempMagMap->d[i*(uint32_t)m_Geometry->N_x+j] = MagUpdateMap->getValue(i, j);
    }

  uint16_t percentile_index=ArrLength/NUM_NON_HOMOGENOUS_ITER;
  //Partial selection sort

  Real_t max;
  uint32_t max_index;
  for(uint32_t i=0; i <= percentile_index;i++)
    {
      max=TempMagMap->d[i];
      max_index=i;
      for(uint32_t j=i+1;j<ArrLength;j++)
        {
    if(TempMagMap->d[j] > max)
      {
                  max=TempMagMap->d[j];
                  max_index=j;
      }
        }
      Real_t temp=TempMagMap->d[i];
      TempMagMap->d[i]=TempMagMap->d[max_index];
      TempMagMap->d[max_index]=temp;
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



//Updates a line of voxels
//Required Global vars / Private members
//ErrorSino
//Weights
//m_Geometry->Object - This is alrady available
//TempCol
//VoxelLine
//NuisanceParams


/* - These get initialized to Zero each time so they are essentially local variables */
//NEIGHBORHOOD
//BOUNDARYFLAG




uint8_t SOCEngine::updateVoxels(int16_t OuterIter, int16_t Iter,
                             VoxelUpdateType updateType,
                             UInt8Image_t::Pointer VisitCount,
                             RNGVars* RandomNumber,
                             AMatrixCol** TempCol,
                             RealVolumeType::Pointer ErrorSino,
                             RealVolumeType::Pointer Weight,
                             AMatrixCol* VoxelLineResponse,
                             ScaleOffsetParams* NuisanceParams,
                             UInt8Image_t::Pointer Mask,
                             CostData::Pointer cost)
{
  uint8_t exit_status=1;//Indicates normal exit ; else indicates to stop inner iterations
  uint16_t subIterations = 1;
  std::string indent("    ");
  uint8_t err = 0;
  uint32_t zero_count = 0;
  Real_t UpdatedVoxelValue;
//  Real_t accuracy=1e-8;
//  uint16_t binarysearch_count
  Real_t accuracy = 1e-9; //This is the rooting accuracy for x
  uint32_t binarysearch_count=10;//Accuracy is 1/(2^10)
  int32_t errorcode = -1;
//  int16_t Idx;

  //FIXME: Where are these Initialized? Or what values should they be initialized to?
  Real_t low = 0.0, high = 0.0;

  if(updateType == RegularRandomOrderUpdate)
  {
    std::cout << indent << "Regular Random Order update of Voxels" << std::endl;
  }
  else if(updateType == HomogeniousUpdate)
  {
    std::cout << indent << "Homogenous update of voxels" << std::endl;
  }
  else if(updateType == NonHomogeniousUpdate)
  {
    std::cout << indent << "Non Homogenous update of voxels" << std::endl;
    subIterations = NUM_NON_HOMOGENOUS_ITER;
  }
  else
  {
    std::cout << indent << "Unknown Voxel Update Type. Returning Now" << std::endl;
    return exit_status;
  }

  Real_t NH_Threshold = 0.0;
  std::stringstream ss;
  int totalLoops = m_TomoInputs->NumOuterIter * m_TomoInputs->NumIter;

  for (uint16_t NH_Iter = 0; NH_Iter < subIterations; ++NH_Iter)
  {
    ss.str("");
    ss << "Outer Iteration: " << OuterIter << " of " << m_TomoInputs->NumOuterIter;
    ss << "   Inner Iteration: " << Iter << " of " << m_TomoInputs->NumIter;
    ss << "   SubLoop: " << NH_Iter << " of " << subIterations;
    float currentLoop = static_cast<float>(OuterIter * m_TomoInputs->NumIter + Iter);
    notify(ss.str(), currentLoop / totalLoops * 100.0f, Observable::UpdateProgressValueAndMessage);
    if(updateType == NonHomogeniousUpdate)
    {
      //Compute VSC and create a map of pixels that are above the threshold value
      ComputeVSC();
      START_TIMER;
      NH_Threshold = SetNonHomThreshold();
      STOP_TIMER;
      PRINT_TIME("  SetNonHomThreshold");
      std::cout << indent <<"NHICD Threshold: "<<NH_Threshold<<std::endl;
      //Use  FiltMagUpdateMap  to find MagnitudeUpdateMask
      //std::cout << "Completed Calculation of filtered magnitude" << std::endl;
      //Calculate the threshold for the top ? % of voxel updates
    }

    //printf("Iter %d\n",Iter);
#ifdef ROI
    //variables used to stop the process
    Real_t AverageUpdate = 0;
    Real_t AverageMagnitudeOfRecon = 0;
#endif

#ifdef RANDOM_ORDER_UPDATES
    int32_t ArraySize = m_Geometry->N_x * m_Geometry->N_z;
    size_t dims[3] =
    { ArraySize, 0, 0 };
    Int32ArrayType::Pointer Counter = Int32ArrayType::New(dims, "Counter");

    for (int32_t j_new = 0; j_new < ArraySize; j_new++)
    {
      Counter->d[j_new] = j_new;
    }
    uint32_t NumVoxelsToUpdate=0;
    for (int32_t j = 0; j < m_Geometry->N_z; j++)
    {
      for (int32_t k = 0; k < m_Geometry->N_x; k++)
      {
        if(updateType == NonHomogeniousUpdate)
        {
          if(MagUpdateMap->getValue(j, k) > NH_Threshold)
          {
            MagUpdateMask->setValue(1, j, k);
            MagUpdateMap->setValue(0, j, k);
            NumVoxelsToUpdate++;
          }
          else
          {
            MagUpdateMask->setValue(0, j, k);
          }
        }
        else if(updateType == HomogeniousUpdate)
        {
          MagUpdateMap->setValue(0, j, k);
          NumVoxelsToUpdate++;
        }
        else if(updateType == RegularRandomOrderUpdate)
        {
          MagUpdateMap->setValue(0, j, k);
          NumVoxelsToUpdate++;
        }
        VisitCount->setValue(0, j, k);
       // VisitCount->d[j][k] = 0;
      }
    }
     std::cout << indent <<"Number of voxel lines to update: "<<NumVoxelsToUpdate<<std::endl;

#endif

    START_TIMER;
    for (int32_t j = 0; j < m_Geometry->N_z; j++) //Row index
    {
      for (int32_t k = 0; k < m_Geometry->N_x; k++) //Column index
      {

#ifdef RANDOM_ORDER_UPDATES
        //RandomNumber=init_genrand(Iter);
        //int32_t Index = (genrand_int31(RandomNumber)) % ArraySize;
        uint32_t Index = (genrand_int31(RandomNumber)) % ArraySize;
//      int32_t k_new = Counter->d[Index] % m_Geometry->N_x;
//      int32_t j_new = Counter->d[Index] / m_Geometry->N_x;
        int32_t k_new = Counter->d[Index] % m_Geometry->N_x;
        int32_t j_new = Counter->d[Index] / m_Geometry->N_x;
        // std::cout<<k_new<<","<<j_new<<std::endl;
        //memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index-1));
        //TODO: Instead just swap the value in Index with the one in ArraySize
        Counter->d[Index] = Counter->d[ArraySize - 1];
//        VisitCount->d[j_new][k_new] = 1;
        VisitCount->setValue(1, j_new, k_new);
        ArraySize--;
        Index = j_new * m_Geometry->N_x + k_new; //This index pulls out the apprppriate index corresponding to
        //the voxel line (j_new,k_new)
#else
        int32_t j_new=j;
        int32_t k_new=k;
        uint32_t Index = j_new*m_Geometry->N_x + k_new;
#endif //Random order updates
        //   AMatrixCol* TempMemBlock = TempCol[j_new][k_new]; //Remove this

        int shouldInitNeighborhood = 0;

        if(updateType == NonHomogeniousUpdate && MagUpdateMask->getValue(j_new, k_new) == 1 && TempCol[Index]->count > 0)
        {
          ++shouldInitNeighborhood;
        }
        if(updateType == HomogeniousUpdate && TempCol[Index]->count > 0)
        {
          ++shouldInitNeighborhood;
        }
        if(updateType == RegularRandomOrderUpdate && TempCol[Index]->count > 0)
        {
          ++shouldInitNeighborhood;
        }

        if(shouldInitNeighborhood > 0)
        //After this should ideally call UpdateVoxelLine(j_new,k_new) ie put everything in this "if" inside a method called UpdateVoxelLine
        {
          for (int32_t i = 0; i < m_Geometry->N_y; i++) //slice index
          {
            //Neighborhood of (i,j,k) should be initialized to zeros each time
            for (int32_t p = 0; p <= 2; p++)
            {
              for (int32_t q = 0; q <= 2; q++)
              {
                for (int32_t r = 0; r <= 2; r++)
                {
                  NEIGHBORHOOD[p][q][r] = 0.0;
                  BOUNDARYFLAG[p][q][r] = 0;
                }
              }
            }
#ifdef CIRCULAR_BOUNDARY_CONDITION
            for(p = -1; p <=1; p++)
            {
              for(q = -1; q <= 1; q++)
              {
                for(r = -1; r <= 1;r++)
                {
                  tempindex_x = mod(r+k_new,m_Geometry->N_x);
                  tempindex_y =mod(p+i,m_Geometry->N_y);
                  tempindex_z = mod(q+j_new,m_Geometry->N_z);
                  NEIGHBORHOOD[p+1][q+1][r+1] = m_Geometry->Object->d[tempindex_z][tempindex_x][tempindex_y];
                  BOUNDARYFLAG[p+1][q+1][r+1]=1;
                }}}
#else
            //For a given (i,j,k) store its 26 point neighborhood
            for (int32_t p = -1; p <= 1; p++)
            {
              for (int32_t q = -1; q <= 1; q++)
              {
                for (int32_t r = -1; r <= 1; r++)
                {
                  if(i + p >= 0 && i + p < m_Geometry->N_y)
                  {
                    if(j_new + q >= 0 && j_new + q < m_Geometry->N_z)
                    {
                      if(k_new + r >= 0 && k_new + r < m_Geometry->N_x)
                      {
                        NEIGHBORHOOD[p + 1][q + 1][r + 1] = m_Geometry->Object->getValue(q + j_new, r + k_new, p + i);
                        BOUNDARYFLAG[p + 1][q + 1][r + 1] = 1;
                      }
                      else
                      {
                        BOUNDARYFLAG[p + 1][q + 1][r + 1] = 0;
                      }
                    }
                  }
                }
              }
            }
#endif//circular boundary condition check
            NEIGHBORHOOD[1][1][1] = 0.0;
#ifndef NDEBUG
            if(i == 0 && j == 31 && k == 31)
            {
              printf("***************************\n");
              printf("Geom %lf\n", m_Geometry->Object->getValue(i, 31, 31) );
              for (int p = 0; p <= 2; p++)
              {
                for (int q = 0; q <= 2; q++)
                {
                  for (int r = 0; r <= 2; r++)
                  {
                    printf("%lf\n", NEIGHBORHOOD[p][q][r]);
                  }
                }
              }
            }
#endif
            //Compute theta1 and theta2
            V = m_Geometry->Object->getValue(j_new, k_new, i); //Store the present value of the voxel
            THETA1 = 0.0;
            THETA2 = 0.0;
#ifdef ZERO_SKIPPING
			  //Zero Skipping Algorithm
			  bool ZSFlag=true;
			  if(V == 0.0 && (Iter > 0 || OuterIter > 0))
			  {
				  for(uint8_t p = 0; p <=2; p++)
					  for(uint8_t q = 0; q <= 2; q++)
						  for(uint8_t r = 0; r <= 2;r++)
							  if(NEIGHBORHOOD[p][q][r] > 0.0)
							  {
								  ZSFlag = false;
								  break;
							  }
			  }
			  else
			  {
				  ZSFlag = false;//First time dont care for zero skipping
			  }
#else
			  bool ZSFlag = false; //do ICD on all voxels
#endif //Zero skipping
			if(ZSFlag == false)
			{
//            Real_t tt;
            //TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
           /* OLD for (uint32_t q = 0; q < TempCol[Index]->count; q++)
            {
              uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (m_Sinogram->N_r)));
              uint16_t i_r = (TempCol[Index]->index[q] % (m_Sinogram->N_r));
              uint16_t VoxelLineAccessCounter = 0;
              Real_t ttmp = NuisanceParams->I_0->d[i_theta]  * TempCol[Index]->values[q];
              for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
              {
                Real_t ProjectionEntry = ttmp * VoxelLineResponse[i].values[VoxelLineAccessCounter];
                tt = Weight->getValue(i_theta, i_r, i_t) * ProjectionEntry;
                THETA2 += (tt * ProjectionEntry);
                THETA1 += (tt * ErrorSino->getValue(i_theta, i_r, i_t) );
                VoxelLineAccessCounter++;
              }
            }*/

			for (uint32_t q = 0; q < TempCol[Index]->count; q++)
			{
					uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (m_Sinogram->N_r)));
					uint16_t i_r = (TempCol[Index]->index[q] % (m_Sinogram->N_r));
					uint16_t VoxelLineAccessCounter = 0;
					for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
					{
                     //size_t weight_idx = Weight->calcIndex(i_theta, i_r, i_t);
                     size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);

					 if (m_Sinogram->BF_Flag == false)
					 {
							Real_t ProjectionEntry = NuisanceParams->I_0->d[i_theta]*VoxelLineResponse[i].values[VoxelLineAccessCounter] * (TempCol[Index]->values[q]);
							THETA2 += (ProjectionEntry*ProjectionEntry*Weight->d[error_idx]);
							THETA1 +=  (ErrorSino->d[error_idx]*ProjectionEntry* Weight->d[error_idx]);
							VoxelLineAccessCounter++;
					  }
					  else
					  {
							//size_t bfcounts_idx = m_BFSinogram->counts->calcIndex(i_theta, i_r, i_t);
							Real_t ProjectionEntry = m_BFSinogram->counts->d[error_idx]*NuisanceParams->I_0->d[i_theta]*VoxelLineResponse[i].values[VoxelLineAccessCounter] * (TempCol[Index]->values[q]);
							THETA2 += (ProjectionEntry*ProjectionEntry*Weight->d[error_idx]);
							THETA1 +=  (ErrorSino->d[error_idx]*ProjectionEntry* Weight->d[error_idx]);
							VoxelLineAccessCounter++;
					  }
					}
				}

            THETA1 *= -1;
            minMax(&low, &high);

#ifdef DEBUG
            if(i == 0 && j == 31 && k == 31) printf("(%lf,%lf,%lf) \n", low, high, V - (THETA1 / THETA2));
#endif

            //Solve the 1-D optimization problem
            //printf("V before updating %lf",V);
#ifndef SURROGATE_FUNCTION
            //TODO : What if theta1 = 0 ? Then this will give error

            DerivOfCostFunc docf(BOUNDARYFLAG, NEIGHBORHOOD, FILTER, V, THETA1, THETA2, SIGMA_X_P, MRF_P);
            UpdatedVoxelValue = (Real_t)solve<DerivOfCostFunc>(&docf, (double)low, (double)high, (double)accuracy, &errorcode,binarysearch_count);

            //std::cout<<low<<","<<high<<","<<UpdatedVoxelValue<<std::endl;
#else
            errorcode = 0;
#ifdef QGGMRF
            UpdatedVoxelValue = CE_FunctionalSubstitution(low, high);
#else
            SurrogateUpdate = surrogateFunctionBasedMin();
            UpdatedVoxelValue = SurrogateUpdate;
#endif //QGGMRF
#endif//Surrogate function
            //printf("%lf\n",SurrogateUpdate);
            if(errorcode == 0)
            {
              //    printf("(%lf,%lf,%lf)\n",low,high,UpdatedVoxelValue);
              //  printf("Updated %lf\n",UpdatedVoxelValue);
#ifdef POSITIVITY_CONSTRAINT
              if(UpdatedVoxelValue < 0.0)
              { //Enforcing positivity constraints
                UpdatedVoxelValue = 0.0;
              }
#endif
            }
            else
            {
              if(THETA1 == 0 && low == 0 && high == 0) UpdatedVoxelValue = 0;
              else
              {
                // printf("Error \n");
                // printf("%d %d\n", j_new, k_new);
              }
            }

            //TODO Print appropriate error messages for other values of error code
            m_Geometry->Object->setValue(UpdatedVoxelValue, j_new, k_new, i);
            //#ifdef NHICD
            Real_t intermediate = MagUpdateMap->getValue(j_new, k_new) + fabs(m_Geometry->Object->getValue(j_new, k_new, i) - V);
            MagUpdateMap->setValue(intermediate, j_new, k_new);
            //#endif

#ifdef ROI
            //if(Mask->d[j_new][k_new] == 1)
            if (Mask->getValue(j_new, k_new) == 1)
            {
              AverageUpdate += fabs(m_Geometry->Object->getValue(j_new, k_new, i) - V);
              AverageMagnitudeOfRecon += fabs(V); //computing the percentage update =(Change in mag/Initial magnitude)
            }
#endif

            //Update the ErrorSinogram -OLD
           /* for (uint32_t q = 0; q < TempCol[Index]->count; q++)
            {
              uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (m_Sinogram->N_r)));
              uint16_t i_r = (TempCol[Index]->index[q] % (m_Sinogram->N_r));
              uint16_t VoxelLineAccessCounter = 0;
              uint32_t index = VoxelLineResponse[i].index[0];
              uint32_t count = VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count;
              for (uint32_t i_t = index; i_t < count; i_t++)
              //for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
              {
                Real_t ttmp = NuisanceParams->I_0->d[i_theta] * (TempCol[Index]->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (m_Geometry->Object->getValue(j_new, k_new, i) - V));
                ErrorSino->deleteFromValue(ttmp, i_theta, i_r, i_t);
                VoxelLineAccessCounter++;
              }
            }*/

				//Update the ErrorSinogram
				for (uint32_t q = 0; q < TempCol[Index]->count; q++)
				{
					uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (m_Sinogram->N_r)));
					uint16_t i_r = (TempCol[Index]->index[q] % (m_Sinogram->N_r));
					uint16_t VoxelLineAccessCounter = 0;
					for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
					{
        		        size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);

						//for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
						if(m_Sinogram->BF_Flag == false)
						{
							ErrorSino->d[error_idx] -= (NuisanceParams->I_0->d[i_theta]
																* (TempCol[Index]->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (m_Geometry->Object->getValue(j_new, k_new, i) - V)));
							VoxelLineAccessCounter++;
						}
						else {
							//size_t bfcounts_idx = m_BFSinogram->counts->calcIndex(i_theta, i_r, i_t);
							ErrorSino->d[error_idx] -= (NuisanceParams->I_0->d[i_theta]*m_BFSinogram->counts->d[error_idx]* (TempCol[Index]->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (m_Geometry->Object->getValue(j_new, k_new, i) - V)));
							VoxelLineAccessCounter++;
						}
					}
				}
		  }
			  else {
				  zero_count++;
			  }

//            Idx++;
          }
        }
        else
        {
          continue;
        }

      }
    }
    STOP_TIMER;
    PRINT_TIME("Voxel Update");

#ifdef RANDOM_ORDER_UPDATES
    for (int j = 0; j < m_Geometry->N_z; j++)
    { //Row index
      for (int k = 0; k < m_Geometry->N_x; k++)
      { //Column index
        //if(VisitCount->d[j][k] == 0)
        if (VisitCount->getValue(j,k) == 0)
        {
          printf("Pixel (%d %d) not visited\n", j, k);
        }
      }
    }
#endif

#ifdef COST_CALCULATE

    /*********************Cost Calculation*************************************/
    Real_t cost_value = computeCost(ErrorSino, Weight);
    std::cout<<cost_value<<std::endl;
    int increase = cost->addCostValue(cost_value);
    if (increase ==1)
    {
      std::cout << "Cost just increased after ICD!" << std::endl;
      break;
    }
    cost->writeCostValue(cost_value);
    /**************************************************************************/
#else
    printf("%d\n",Iter);
#endif //Cost calculation endif
#ifdef ROI
    if(AverageMagnitudeOfRecon > 0)
    {
      printf("%d,%lf\n", Iter + 1, AverageUpdate / AverageMagnitudeOfRecon);

    //Use the stopping criteria if we are performing a full update of all voxels
      if((AverageUpdate / AverageMagnitudeOfRecon) < m_TomoInputs->StopThreshold && updateType != NonHomogeniousUpdate)
      {
        printf("This is the terminating point %d\n", Iter);
        m_TomoInputs->StopThreshold *= THRESHOLD_REDUCTION_FACTOR; //Reducing the thresold for subsequent iterations
    std::cout<<"New threshold"<<m_TomoInputs->StopThreshold<<std::endl;
      exit_status=0;
        break;
      }
    }
#endif//ROI end

#ifdef WRITE_INTERMEDIATE_RESULTS

    if(Iter == NumOfWrites*WriteCount)
    {
      WriteCount++;
      sprintf(buffer,"%d",Iter);
      sprintf(Filename,"ReconstructedObjectAfterIter");
      strcat(Filename,buffer);
      strcat(Filename,".bin");
      Fp3 = fopen(Filename, "w");
      //  for (i=0; i < Geometry->N_y; i++)
      //    for (j=0; j < Geometry->N_z; j++)
      //      for (k=0; k < Geometry->N_x; k++)
      TempPointer = m_Geometry->Object;
      NumOfBytesWritten=fwrite(&(m_Geometry->Object->d[0][0][0]), sizeof(Real_t),m_Geometry->N_x*m_Geometry->N_y*m_Geometry->N_z, Fp3);
      printf("%d\n",NumOfBytesWritten);

      fclose(Fp3);
    }
#endif

    if(getCancel() == true)
    {
      setErrorCondition(err);
      return exit_status;
    }

  }

	std::cout<<"Number of zeroed out entries="<<zero_count<<std::endl;
  return exit_status;

}


#ifdef QGGMRF
//Function to compute parameters of thesurrogate function
void SOCEngine::CE_ComputeQGGMRFParameters(Real_t umin,Real_t umax,Real_t RefValue)
{
//  DATA_TYPE DeltaMin,DeltaMax,T,Derivative_Delta0;
   Real_t Delta0;
//  DATA_TYPE AbsDelta0,AbsDeltaMin,AbsDeltaMax;
  uint8_t i,j,k,count=0;
  for(i=0;i<3;i++)
    for (j=0; j < 3; j++)
      for(k=0; k < 3;k++)
        if ((i != 1 || j !=1 || k != 1) &&  BOUNDARYFLAG[i][j][k] == 1)
      /*  {
          Delta0 = V - NEIGHBORHOOD[i][j][k];
          DeltaMin = umin - NEIGHBORHOOD[i][j][k];
          DeltaMax = umax - NEIGHBORHOOD[i][j][k];
          AbsDelta0 = fabs(Delta0);
          AbsDeltaMin = fabs(DeltaMin);
          AbsDeltaMax = fabs(DeltaMax);

          if(AbsDelta0 <= Minimum(AbsDeltaMin,AbsDeltaMax))
            T = -Delta0;
          else if(AbsDeltaMin <= Minimum(AbsDelta0,AbsDeltaMax))
            T = DeltaMin;
          else if(AbsDeltaMax <= Minimum(AbsDeltaMin,AbsDelta0))
            T = DeltaMax;
          if(Delta0 != 0)
          {
            Derivative_Delta0 = CE_QGGMRF_Derivative(Delta0);
            if( T == -Delta0)
            {
            QGGMRF_Params[count][0] = - Derivative_Delta0/(T-Delta0); //Because the QGGMRF prior is symmetric
            }
            else
            {
            QGGMRF_Params[count][0] = ((CE_QGGMRF_Value(T) - CE_QGGMRF_Value(Delta0))/((T-Delta0)*(T - Delta0))) - (Derivative_Delta0/(T-Delta0));
            }

          }
          else
          {
            Derivative_Delta0 = CE_QGGMRF_Derivative(Delta0);
            QGGMRF_Params[count][0] = CE_QGGMRF_SecondDerivative(0)/2;
          }
          QGGMRF_Params[count][1] = Derivative_Delta0 - 2*QGGMRF_Params[count][0]*V;
          QGGMRF_Params[count][2] = CE_QGGMRF_Value(Delta0) - QGGMRF_Params[count][0]*V*V - QGGMRF_Params[count][1]*V;
          count++;
        }*/
        {
          Delta0  = RefValue - NEIGHBORHOOD[i][j][k];

          if(Delta0 != 0)
          QGGMRF_Params[count][0] = CE_QGGMRF_Derivative(Delta0)/(Delta0);
          else {
            QGGMRF_Params[count][0] = CE_QGGMRF_SecondDerivative(0);
          }
          //QGGMRF_Params[count][1] = CE_QGGMRF_Value(Delta0) - (Delta0/2)*CE_QGGMRF_Derivative(Delta0);
          count++;
        }

}

Real_t SOCEngine::CE_FunctionalSubstitution(Real_t umin,Real_t umax)
{
  Real_t u,temp1=0,temp2=0,temp_const,RefValue=0;
  uint8_t i,j,k,count=0;
#ifdef POSITIVITY_CONSTRAINT
  if(umin < 0)
    umin =0;
#endif //Positivity
  RefValue = V;
  //Need to Loop this for multiple iterations of substitute function
  for(int8_t Iter=0; Iter < QGGMRF_ITER;Iter++)
  {
  CE_ComputeQGGMRFParameters(umin, umax,RefValue);
  for(i = 0;i < 3; i++)
    for (j=0; j < 3; j++)
      for(k=0; k < 3; k++)
      {
        if((i != 1 || j != 1 || k != 1) && BOUNDARYFLAG[i][j][k] == 1)
        /*{
          temp1 += FILTER[i][j][k]*QGGMRF_Params[count][1];
          temp2 += FILTER[i][j][k]*QGGMRF_Params[count][0];
            count++;
        }*/
        {
          temp_const = FILTER[i][j][k]*QGGMRF_Params[count][0];
          temp1 += temp_const*NEIGHBORHOOD[i][j][k];
          temp2 += temp_const;
          count++;
        }
      }
  u=(temp1+ (THETA2*V) - THETA1)/(temp2 + THETA2);

    if(Iter < QGGMRF_ITER-1)
      RefValue = Clip(RefValue + MRF_ALPHA*(u-RefValue),umin,umax);
    else {
		RefValue = Clip(u, umin, umax);
    }

  }

  return RefValue;

}




Real_t SOCEngine::CE_QGGMRF_Value(Real_t delta)
{
  return ((pow(fabs(delta),MRF_P)/SIGMA_X_P)/(MRF_C + pow(fabs(delta),MRF_P - MRF_Q)/SIGMA_X_P_Q));
}

Real_t SOCEngine::CE_QGGMRF_Derivative(Real_t delta)
{
  Real_t temp1,temp2,temp3;
  temp1=pow(fabs(delta),MRF_P - MRF_Q)/(SIGMA_X_P_Q);
  temp2=pow(fabs(delta),MRF_P - 1);
  temp3 = MRF_C + temp1;
  if(delta < 0)
    return ((-1*temp2/(temp3*SIGMA_X_P))*(MRF_P - ((MRF_P-MRF_Q)*temp1)/(temp3)));
  else
  {
    return ((temp2/(temp3*SIGMA_X_P))*(MRF_P - ((MRF_P-MRF_Q)*temp1)/(temp3)));
  }
}

Real_t SOCEngine::CE_QGGMRF_SecondDerivative(Real_t delta)
{
  return MRF_P/(SIGMA_X_P*MRF_C);
}

#endif //QGGMRF

#include "SOCEngine_Extra.cpp"


