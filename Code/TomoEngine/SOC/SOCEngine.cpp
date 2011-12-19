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
#include "TomoEngine/IO/RawGeometryWriter.h"
#include "TomoEngine/IO/NuisanceParamWriter.h"
#include "TomoEngine/IO/SinogramBinWriter.h"
#include "TomoEngine/IO/VTKFileWriters.hpp"
#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/SOC/ForwardProject.h"
#include "TomoEngine/Filters/ComputeGainsOffsets.h"
#include "TomoEngine/Filters/DetectorResponse.h"
#include "TomoEngine/Filters/DetectorResponseWriter.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/Filters/MRCSinogramInitializer.h"
#include "TomoEngine/Filters/RawSinogramInitializer.h"
#include "TomoEngine/Filters/GainsOffsetsReader.h"
#include "TomoEngine/Filters/ComputeGainsOffsets.h"
#include "TomoEngine/Filters/InitialReconstructionInitializer.h"
#include "TomoEngine/Filters/InitialReconstructionBinReader.h"




#define START_TIMER startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME std::cout << indent << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

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
SOCEngine::SOCEngine() :
    m_Inputs(NULL),
    m_Sinogram(NULL),
    m_Geometry(NULL),
    m_NuisanceParams(NULL)
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


  m_Inputs = NULL;
  m_Sinogram = NULL;
  m_Geometry = NULL;
  m_NuisanceParams = NULL;

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
  HAMMING_WINDOW[2][0]=0.0159;HAMMING_WINDOW[2][1]=0.1076;HAMMING_WINDOW[2][2]=0.1993;HAMMING_WINDOW[2][3]=0.1076;HAMMING_WINDOW[2][4]=0.0159;
  HAMMING_WINDOW[3][0]= 0.0013;HAMMING_WINDOW[3][1]=0.0086;HAMMING_WINDOW[3][2]=0.0159;HAMMING_WINDOW[3][3]=0.0086;HAMMING_WINDOW[3][4]=0.0013;
  HAMMING_WINDOW[4][0]= 0.0086;HAMMING_WINDOW[4][1]=0.0581;HAMMING_WINDOW[4][2]=0.1076;HAMMING_WINDOW[4][3]=0.0581;HAMMING_WINDOW[4][4]=0.0086;


}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::execute()
{
   uint64_t totalTime = EIMTOMO_getMilliSeconds();


	int16_t i,j,k,Idx;
	size_t dims[3];
	int16_t index_delta_t;

  Int32ArrayType::Pointer Counter;
  UInt8ImageType::Pointer VisitCount;

	//uint16_t cost_counter=0;

	uint16_t MaxNumberOfDetectorElts;

	DATA_TYPE center_t,delta_t;
	DATA_TYPE w3,w4;
	DATA_TYPE checksum = 0,temp;

	RealImageType::Pointer VoxelProfile;
	RealVolumeType::Pointer detectorResponse;
	RealVolumeType::Pointer H_t;
	DATA_TYPE ProfileCenterT;


	RealVolumeType::Pointer Y_Est;//Estimated Sinogram
	RealVolumeType::Pointer Final_Sinogram;//To store and write the final sinogram resulting from our reconstruction
	RealVolumeType::Pointer ErrorSino;//Error Sinogram
	RealVolumeType::Pointer Weight;//This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation
	RNGVars* RandomNumber;
	std::string indent("");

#if TomoEngine_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
#endif

  // Initialize the Sinogram
  if (m_Inputs == NULL)
  {
    setErrorCondition(-1);
    updateProgressAndMessage("Error: The TomoInput Structure was NULL. The proper API is to supply this class with that structure,", 100);
    return;
  }
  //Based on the inputs , calculate the "other" variables in the structure definition
  if (m_Sinogram == NULL)
  {
    setErrorCondition(-1);
    updateProgressAndMessage("Error: The Sinogram Structure was NULL. The proper API is to supply this class with that structure,", 100);
    return;
  }
  if (m_Geometry == NULL)
  {
    setErrorCondition(-1);
    updateProgressAndMessage("Error: The Geometry Structure was NULL. The proper API is to supply this class with that structure,", 100);
    return;
  }



	// Read the Input data from the supplied data file
	// We are scoping here so the various readers are automatically cleaned up before
	// the code goes any farther
	{
	  TomoFilter::Pointer dataReader = TomoFilter::NullPointer();
    std::string extension = MXAFileInfo::extension(m_Inputs->SinoFile);
    if (extension.compare("mrc") == 0 || extension.compare("ali") == 0)
    {
      dataReader = MRCSinogramInitializer::NewTomoFilter();
    }
    else if (extension.compare("bin") == 0 )
    {
      dataReader = RawSinogramInitializer::NewTomoFilter();
    }
    else
    {
      setErrorCondition(-1);
      updateProgressAndMessage("A supported file reader for the input file was not found." , 100);
      return;
    }
    dataReader->setInputs(m_Inputs);
    dataReader->setSinogram(m_Sinogram);
    dataReader->addObserver(this);
    dataReader->execute();
    if (dataReader->getErrorCondition() < 0)
    {
      updateProgressAndMessage("Error reading Input Sinogram Data file", 100);
      setErrorCondition(dataReader->getErrorCondition());
      return;
    }
	}


  // Now read or generate the Gains and Offsets data. We are scoping this section
  // so the reader automactically gets cleaned up at this point.
  {
    TomoFilter::Pointer gainsOffsetsInitializer = TomoFilter::NullPointer();
    if(m_Inputs->GainsOffsetsFile.empty() == true)      // Calculate the initial Gains and Offsets
    {
      gainsOffsetsInitializer = ComputeGainsOffsets::NewTomoFilter();
    }
    else
    {
      gainsOffsetsInitializer = GainsOffsetsReader::NewTomoFilter();
    }
    gainsOffsetsInitializer->setSinogram(m_Sinogram);
    gainsOffsetsInitializer->setInputs(m_Inputs);
    gainsOffsetsInitializer->addObserver(this);
    gainsOffsetsInitializer->execute();
    if(gainsOffsetsInitializer->getErrorCondition() < 0)
    {
     updateProgressAndMessage("Error initializing Input Gains and Offsets Data file", 100);
     setErrorCondition(gainsOffsetsInitializer->getErrorCondition());
     return;
    }
    /********************REMOVE************************/
    std::cout<<"HARD WIRED TARGET GAIN"<<std::endl;

    m_Sinogram->TargetGain = TARGET_GAIN;
    std::cout << "Target Gain: " << m_Sinogram->TargetGain << std::endl;
    /*************************************************/
  }

  // Initialize the Geometry data from a rough reconstruction
  {
    InitialReconstructionInitializer::Pointer geomInitializer = InitialReconstructionInitializer::NullPointer();
    std::string extension = MXAFileInfo::extension(m_Inputs->InitialReconFile);
    if (m_Inputs->InitialReconFile.empty() == true)
    {
      // This will just initialize all the values to Zero (0)
      geomInitializer = InitialReconstructionInitializer::New();
    }
    else if (extension.compare("bin") == 0 )
    {
      // This will read the values from a binary file
      geomInitializer = InitialReconstructionBinReader::NewInitialReconstructionInitializer();
    }
    geomInitializer->setSinogram(m_Sinogram);
    geomInitializer->setInputs(m_Inputs);
    geomInitializer->setGeometry(m_Geometry);
    geomInitializer->addObserver(this);
    geomInitializer->execute();

    if(geomInitializer->getErrorCondition() < 0)
    {
      updateProgressAndMessage("Error reading Initial Reconstruction Data from File", 100);
      setErrorCondition(geomInitializer->getErrorCondition());
      return;
    }
  }

	//Scale and Offset Parameter Structures
	ScaleOffsetParams NuisanceParams;
	NuisanceParams.I_0 = RealArrayType::NullPointer();
	NuisanceParams.mu = RealArrayType::NullPointer();
	NuisanceParams.alpha = RealArrayType::NullPointer();

//	DATA_TYPE numerator_sum =0;
//	DATA_TYPE denominator_sum = 0;
	DATA_TYPE x,z;
	DATA_TYPE sum;
//	DATA_TYPE a,b,c,d,e;
//	DATA_TYPE determinant;
//	DATA_TYPE LagrangeMultiplier;
  //This is used to store the projection of the object for each view
  dims[1] = m_Sinogram->N_t;
  dims[0] = m_Sinogram->N_r;
  RealImageType::Pointer MicroscopeImage = RealImageType::New(dims);
  MicroscopeImage->setName("MicroscopeImage");


  dims[0] = m_Sinogram->N_theta;
  NuisanceParams.I_0 = RealArrayType::New(dims);
  NuisanceParams.I_0->setName("NuisanceParams.I_0");
  NuisanceParams.mu = RealArrayType::New(dims);
  NuisanceParams.mu->setName("NuisanceParams.mu");
#ifdef NOISE_MODEL
	//alpha is the noise variance adjustment factor
	NuisanceParams.alpha = RealArrayType::New(dims);
	NuisanceParams.alpha->setName("NuisanceParams.alpha");
#endif

  // initialize variables
  Idx = 0;

#ifdef WRITE_INTERMEDIATE_RESULTS
	DATA_TYPE Fraction = 0.1;//write this fraction of the iterations
	int16_t NumOfWrites = floor((DATA_TYPE)(m_Inputs->NumIter)*Fraction);
	int16_t WriteCount = 0;
	char Filename[100];
	char buffer[20];
	void* TempPointer;
	size_t NumOfBytesWritten;
#endif

#ifdef ROI
	UInt8ImageType::Pointer Mask;
	DATA_TYPE EllipseA,EllipseB;
#endif

#ifdef COST_CALCULATE
  //cost=(DATA_TYPE*)get_spc((m_Inputs->NumIter+1)*m_Inputs->NumOuterIter*3,sizeof(DATA_TYPE));//the factor 3 is in there to ensure we can store 3 costs per iteration one after update x, then mu and then I
  dims[0] = (m_Inputs->NumIter+1)*m_Inputs->NumOuterIter*3;
//  cost = RealArrayType::New(dims);
//  cost->setName("cost");
#endif

	dims[0] = m_Sinogram->N_theta;
	dims[1] = m_Sinogram->N_r;
	dims[2] = m_Sinogram->N_t;

  //  Y_Est=(RealVolumeType::Pointer )get_3D(m_Sinogram->N_theta,m_Sinogram->N_r,m_Sinogram->N_t,sizeof(DATA_TYPE));
  //  ErrorSino=(RealVolumeType::Pointer )get_3D(m_Sinogram->N_theta,m_Sinogram->N_r,m_Sinogram->N_t,sizeof(DATA_TYPE));
  //  Weight=(RealVolumeType::Pointer )get_3D(m_Sinogram->N_theta,m_Sinogram->N_r,m_Sinogram->N_t,sizeof(DATA_TYPE));
  //
  Y_Est = RealVolumeType::New(dims);
  Y_Est->setName("Y_Est");
  ErrorSino = RealVolumeType::New(dims);
  ErrorSino->setName("ErrorSino");
  Weight = RealVolumeType::New(dims);
  Weight->setName("Weight");

  //Setting the value of all the global variables
  OffsetR = ((m_Inputs->delta_xz/sqrt(3.0)) + m_Sinogram->delta_r/2)/DETECTOR_RESPONSE_BINS;
  OffsetT = ((m_Inputs->delta_xz/2) + m_Sinogram->delta_t/2)/DETECTOR_RESPONSE_BINS;

#ifdef BEAM_CALCULATION
	BEAM_WIDTH = (0.5)*m_Sinogram->delta_r;
#else
	BEAM_WIDTH = m_Sinogram->delta_r;
#endif

#ifndef QGGMRF
	MRF_P = m_Inputs->p;
	SIGMA_X_P = pow(m_Inputs->SigmaX,MRF_P);
#else
	MRF_P=2;
	MRF_Q=1.2;
	MRF_C=30;
	MRF_ALPHA=1.5;
	SIGMA_X_P = pow(m_Inputs->SigmaX,MRF_P);
	for(i=0;i < 3;i++) {
		for(j=0; j < 3;j++) {
			for(k=0; k < 3;k++) {
				FILTER[i][j][k]*=(1.0/SIGMA_X_P); }}}
#endif //QGGMRF

	//globals assosiated with finding the optimal gain and offset parameters

	dims[0] = m_Sinogram->N_theta;
	dims[1] = 3;
	dims[2] = 0;
	//Hold the coefficients of a quadratic equation
	//  QuadraticParameters = (DATA_TYPE**)(get_img(3, m_Sinogram->N_theta, sizeof(DATA_TYPE)));
	QuadraticParameters = RealImageType::New(dims);
	QuadraticParameters->setName("QuadraticParameters");
	//  Qk_cost=(DATA_TYPE**)get_img(3, m_Sinogram->N_theta, sizeof(DATA_TYPE));
	Qk_cost = RealImageType::New(dims);
  Qk_cost->setName("Qk_cost");
//	bk_cost=(DATA_TYPE**)get_img(2, m_Sinogram->N_theta, sizeof(DATA_TYPE));
    dims[1] = 2;
	bk_cost = RealImageType::New(dims);
	bk_cost->setName("bk_cost");

	dims[0] = m_Sinogram->N_theta;
	ck_cost = RealArrayType::New(dims);
	ck_cost->setName("ck_cost");
//	ck_cost=(DATA_TYPE*)get_spc(m_Sinogram->N_theta, sizeof(DATA_TYPE));

	NumOfViews = m_Sinogram->N_theta;
	LogGain = m_Sinogram->N_theta*log(m_Sinogram->TargetGain);
	dims[0] = m_Sinogram->N_theta;
//	d1=(DATA_TYPE*)get_spc(m_Sinogram->N_theta, sizeof(DATA_TYPE));
	d1 = RealArrayType::New(dims);
	d1->setName("d1");
//	d2=(DATA_TYPE*)get_spc(m_Sinogram->N_theta, sizeof(DATA_TYPE));
	d2 = RealArrayType::New(dims);
	d2->setName("d2");

	//calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram Structure are converted to radians
	VoxelProfile = calculateVoxelProfile(); //Verified with ML
	//Pre compute sine and cos theta to speed up computations
	calculateSinCos();
	//Initialize the e-beam
	initializeBeamProfile(); //verified with ML

	//calculate sine and cosine of all angles and store in the global arrays sine and cosine
	DetectorResponse::Pointer dResponseFilter = DetectorResponse::New();
	dResponseFilter->setInputs(m_Inputs);
	dResponseFilter->setSinogram(m_Sinogram);
	dResponseFilter->setBeamWidth(BEAM_WIDTH);
	dResponseFilter->setOffsetR(OffsetR);
	dResponseFilter->setOffsetT(OffsetT);
	dResponseFilter->setVoxelProfile(VoxelProfile);
	dResponseFilter->setBeamProfile(BeamProfile);
	dResponseFilter->addObserver(this);
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
	responseWriter->setInputs(m_Inputs);
	responseWriter->setSinogram(m_Sinogram);
	responseWriter->addObserver(this);
	responseWriter->setResponse(detectorResponse);
	responseWriter->execute();
  if (responseWriter->getErrorCondition() < 0)
  {
    std::cout << "Error writing detector response to file." << __FILE__ << "(" << __LINE__ << ")" << std::endl;
    setErrorCondition(-2);
    updateProgressAndMessage("Error Encountered During Reconstruction", 100);
    return;
  }

	dims[0] = 1;
	dims[1] = m_Sinogram->N_theta;
	dims[2] = DETECTOR_RESPONSE_BINS;

	H_t = RealVolumeType::New(dims);
	H_t->setName("H_t");

#ifdef RANDOM_ORDER_UPDATES

  //VisitCount=(uint8_t **)get_img(m_Geometry->N_x, m_Geometry->N_z,sizeof(uint8_t));//width,height
  dims[0] = m_Geometry->N_z;
  dims[1] = m_Geometry->N_x;
  VisitCount = UInt8ImageType::New(dims);
  VisitCount->setName("VisitCount");
// Initialize the Array to zero
  ::memset( VisitCount->d[0], 0, dims[0] * dims[1] * sizeof(uint8_t));
//  for (i = 0; i < m_Geometry->N_z; i++)
//  {
//    for (j = 0; j < m_Geometry->N_x; j++)
//    {
//      VisitCount->d[i][j] = 0;
//    }
//  }
#endif//Random update


  //#ifdef NHICD
	dims[0]=m_Geometry->N_z;//height
	dims[1]=m_Geometry->N_x;//width
	dims[2]=0;

	MagUpdateMap = RealImageType::New(dims);
	MagUpdateMap->setName("Update Map for voxel lines");

	FiltMagUpdateMap = RealImageType::New(dims);
	FiltMagUpdateMap->setName("Update Map for voxel lines");

	MagUpdateMask = Uint8ImageType::New(dims);
	MagUpdateMask->setName("Update Mask for selecting voxel lines NHICD");
	//#endif


#ifdef ROI
  //Mask = (uint8_t**)get_img(m_Geometry->N_x, m_Geometry->N_z,sizeof(uint8_t));//width,height
  dims[0] = m_Geometry->N_z;
  dims[1] = m_Geometry->N_x;
  Mask = UInt8ImageType::New(dims);
  Mask->setName("Mask");
  EllipseA = m_Geometry->LengthX/2;
  EllipseB = m_Inputs->LengthZ/2;
  for (i = 0; i < m_Geometry->N_z ; i++)
  {
    for (j = 0; j < m_Geometry->N_x; j++)
    {
      x = m_Geometry->x0 + ((DATA_TYPE)j + 0.5)*m_Inputs->delta_xz;
      z = m_Geometry->z0 + ((DATA_TYPE)i + 0.5)*m_Inputs->delta_xz;
      if (x >= -(m_Sinogram->N_r*m_Sinogram->delta_r)/2
          && x <= (m_Sinogram->N_r*m_Sinogram->delta_r)/2
          && z >= -m_Inputs->LengthZ/2
          && z <= m_Inputs->LengthZ/2)
      {
        Mask->d[i][j] = 1;
      }
      else
      {
        Mask->d[i][j] = 0;
      }
    }
  }
#endif
	//m_Sinogram->TargetGain=20000;

#ifdef BRIGHT_FIELD //Take log of the data and subtract log(Dosage) from it
	for (i_theta = 0;i_theta < m_Sinogram->N_theta; i_theta++) {//slice index
		for(i_r = 0; i_r < m_Sinogram->N_r;i_r++) {
			for(i_t = 0;i_t < m_Sinogram->N_t;i_t++)
			{

				m_Sinogram->counts->d[i_theta][i_r][i_t] = log(m_Sinogram->counts->d[i_theta][i_r][i_t])-log(m_Sinogram->InitialGain[i_theta]);
			}}}

#endif//Bright Field

  //Scale and Offset Parameters Initialization
  sum = 0;
  temp = 0;
  for(k=0 ; k < m_Sinogram->N_theta;k++)
  {
    // Gains
    NuisanceParams.I_0->d[k] = m_Sinogram->InitialGain->d[k];
    // Offsets
    NuisanceParams.mu->d[k] = m_Sinogram->InitialOffset->d[k];
#ifdef GEOMETRIC_MEAN_CONSTRAINT
    sum+=log(NuisanceParams.I_0->d[k]);
#else
    sum+=NuisanceParams.I_0->d[k];
#endif
  }
	sum/=m_Sinogram->N_theta;
#ifdef GEOMETRIC_MEAN_CONSTRAINT
	sum=exp(sum);
	printf("The geometric mean of the gains is %lf\n",sum);

	//Checking if the input parameters satisfy the target geometric mean
	if(fabs(sum - m_Sinogram->TargetGain) > 1e-5)
	{
		printf("The input paramters dont meet the constraint..Renormalizing\n");
    temp = exp(log(m_Sinogram->TargetGain) - log(sum));
		for (k = 0 ; k < m_Sinogram->N_theta; k++) {
			NuisanceParams.I_0->d[k]=m_Sinogram->InitialGain->d[k]*temp;
		}
	}
#else
	printf("The Arithmetic mean of the constraint is %lf\n",sum);
	if(sum - m_Sinogram->TargetGain > 1e-5)
	{
		printf("Arithmetic Mean Constraint not met..renormalizing\n");
		temp = m_Sinogram->TargetGain/sum;
		for (k = 0 ; k < m_Sinogram->N_theta; k++) {
			NuisanceParams.I_0->d[k]=m_Sinogram->InitialGain->d[k]*temp;
		}
	}
#endif


  for(k = 0 ; k <m_Sinogram->N_theta; k++)
  {
    for (i = 0; i < DETECTOR_RESPONSE_BINS; i++)
    {
      ProfileCenterT = i * OffsetT;
      if(m_Inputs->delta_xy >= m_Sinogram->delta_t)
      {
        if(ProfileCenterT <= ((m_Inputs->delta_xy / 2) - (m_Sinogram->delta_t / 2))) H_t->d[0][k][i] = m_Sinogram->delta_t;
        else
        {
          H_t->d[0][k][i] = -1 * ProfileCenterT + (m_Inputs->delta_xy / 2) + m_Sinogram->delta_t / 2;
        }
        if(H_t->d[0][k][i] < 0)
        {
          H_t->d[0][k][i] = 0;
        }

      }
      else
      {
        if(ProfileCenterT <= m_Sinogram->delta_t / 2 - m_Inputs->delta_xy / 2)
        {
          H_t->d[0][k][i] = m_Inputs->delta_xy;
        }
        else
        {
          H_t->d[0][k][i] = -ProfileCenterT + (m_Inputs->delta_xy / 2) + m_Sinogram->delta_t / 2;
        }

        if(H_t->d[0][k][i] < 0)
        {
          H_t->d[0][k][i] = 0;
        }

      }
    }
  }

	/*for(i=0;i<Sinogram->N_theta;i++)
		for(j=0;j< PROFILE_RESOLUTION;j++)
			checksum+=VoxelProfile[i][j];
    printf("CHK SUM%lf\n",checksum);
*/
    checksum=0;

	//Allocate space for storing columns the A-matrix; an array of pointers to columns
	//AMatrixCol** AMatrix=(AMatrixCol **)get_spc(Geometry->N_x*Geometry->N_z,sizeof(AMatrixCol*));
//	DetectorResponse = CE_DetectorResponse(0,0,Sinogram,Geometry,VoxelProfile);//System response

#ifdef STORE_A_MATRIX

	AMatrixCol**** AMatrix = (AMatrixCol ****)multialloc(sizeof(AMatrixCol*),3,m_Geometry->N_y,m_Geometry->N_z,m_Geometry->N_x);
#else

	//TODO: All this needs to be deallocated at some point
	DATA_TYPE y;
	DATA_TYPE t,tmin,tmax,ProfileThickness;
	int16_t slice_index_min,slice_index_max;
	AMatrixCol*** TempCol = (AMatrixCol***)multialloc(sizeof(AMatrixCol*),2,m_Geometry->N_z,m_Geometry->N_x);//stores 2-D forward projector
//	AMatrixCol* Aj;
//	AMatrixCol* TempMemBlock;
	//T-direction response
	AMatrixCol* VoxelLineResponse=(AMatrixCol*)get_spc(m_Geometry->N_y, sizeof(AMatrixCol));
	MaxNumberOfDetectorElts = (uint16_t)((m_Inputs->delta_xy/m_Sinogram->delta_t)+2);
  for (i=0; i < m_Geometry->N_y; i++) {
		VoxelLineResponse[i].count=0;
		VoxelLineResponse[i].values=(DATA_TYPE*)get_spc(MaxNumberOfDetectorElts, sizeof(DATA_TYPE));
		VoxelLineResponse[i].index=(uint32_t*)get_spc(MaxNumberOfDetectorElts, sizeof(uint32_t));
	}
#endif

	//Calculating A-Matrix one column at a time
	//For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
	// And then store only the non zero entries by allocating a new array of the desired size
	//k=0;

	checksum = 0;
	//q = 0;

#ifdef STORE_A_MATRIX
	for(i = 0;i < m_Geometry->N_y; i++){
		for(j = 0;j < m_Geometry->N_z; j++){
			for (k = 0; k < m_Geometry->N_x; k++)
			{
				//  AMatrix[q++]=CE_CalculateAMatrixColumn(i,j,Sinogram,Geometry,VoxelProfile);
				AMatrix[i][j][k]=calculateAMatrixColumn(j,k,i,m_Sinogram,m_Geometry,VoxelProfile);//row,col,slice

				for(p = 0; p < AMatrix[i][j][k]->count; p++) {
					checksum += AMatrix[i][j][k]->values[p];}
				//   printf("(%d,%d,%d) %lf \n",i,j,k,AMatrix[i][j][k]->values);
				checksum = 0;
			}}}
	printf("Stored A matrix\n");
#else
	temp = 0;
  for (j = 0; j < m_Geometry->N_z; j++)
  {
    for (k = 0; k < m_Geometry->N_x; k++)
    {
      TempCol[j][k] = (AMatrixCol*)calculateAMatrixColumnPartial(j, k, 0, detectorResponse);
      temp += TempCol[j][k]->count;
    }
  }
#endif

	//Storing the response along t-direction for each voxel line
  updateProgressAndMessage("Storing the response along Y-direction for each voxel line", 0);
	for (i =0; i < m_Geometry->N_y; i++)
	{
    y = ((DATA_TYPE)i + 0.5) * m_Inputs->delta_xy + m_Geometry->y0;
    t = y;
    tmin = (t - m_Inputs->delta_xy / 2) > m_Sinogram->T0 ? t - m_Inputs->delta_xy / 2 : m_Sinogram->T0;
    tmax = (t + m_Inputs->delta_xy / 2) <= m_Sinogram->TMax ? t + m_Inputs->delta_xy / 2 : m_Sinogram->TMax;

    slice_index_min = floor((tmin - m_Sinogram->T0) / m_Sinogram->delta_t);
    slice_index_max = floor((tmax - m_Sinogram->T0) / m_Sinogram->delta_t);

    if(slice_index_min < 0)
    {
      slice_index_min = 0;
    }
    if(slice_index_max >= m_Sinogram->N_t)
    {
      slice_index_max = m_Sinogram->N_t - 1;
    }

    //printf("%d %d\n",slice_index_min,slice_index_max);

    for (int i_t = slice_index_min; i_t <= slice_index_max; i_t++)
    {
      center_t = ((DATA_TYPE)i_t + 0.5) * m_Sinogram->delta_t + m_Sinogram->T0;
      delta_t = fabs(center_t - t);
      index_delta_t = floor(delta_t / OffsetT);
      if(index_delta_t < DETECTOR_RESPONSE_BINS)
      {
        w3 = delta_t - (DATA_TYPE)(index_delta_t) * OffsetT;
        w4 = ((DATA_TYPE)index_delta_t + 1) * OffsetT - delta_t;
        ProfileThickness = (w4 / OffsetT) * H_t->d[0][0][index_delta_t]
            + (w3 / OffsetT) * H_t->d[0][0][index_delta_t + 1 < DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : DETECTOR_RESPONSE_BINS - 1];
      }
      else
      {
        ProfileThickness = 0;
      }

      if(ProfileThickness != 0) //Store the response of this slice
      {
        printf("%d %lf\n", i_t, ProfileThickness);
        VoxelLineResponse[i].values[VoxelLineResponse[i].count] = ProfileThickness;
        VoxelLineResponse[i].index[VoxelLineResponse[i].count++] = i_t;
      }
    }
  }


	printf("Number of non zero entries of the forward projector is %lf\n",temp);
	printf("Geometry-Z %d\n",m_Geometry->N_z);


	//Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
	//is Y_Est initailized to zero?
	for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (j = 0; j < m_Sinogram->N_r; j++)
    {
      for (k = 0; k < m_Sinogram->N_t; k++)
      {
        Y_Est->d[i][j][k] = 0.0;
      }
    }
  }


	RandomNumber=init_genrand(1ul);

  updateProgressAndMessage("Starting Forward Projection", 10);
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
    ForwardProject fp(m_Sinogram, m_Geometry, TempCol,VoxelLineResponse,Y_Est, &NuisanceParams, t);
    fp();
#endif

  }
#if TomoEngine_USE_PARALLEL_ALGORITHMS
  g->wait(); // Wait for all the threads to complete before moving on.
  delete g;
#endif
  printf("\n");
  STOP_TIMER;
  printf("Forward Project Time was ");
  PRINT_TIME;


  //Calculate Error Sinogram - Can this be combined with previous loop?
  //Also compute weights of the diagonal covariance matrix
#ifdef NOISE_MODEL
  START_TIMER;

  //TODO: This can be parallelized  for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++) //slice index
  for (int16_t zz = 0; zz < m_Sinogram->N_theta; zz++) //slice index
  {
    std::cout << "\rComputing Weights for Tilt Index " << zz<< "/" << m_Sinogram->N_theta;
    NuisanceParams.alpha->d[zz] = 1; //Initialize the refinement parameters to 1
    for (uint16_t xx = 0; xx < m_Sinogram->N_r; xx++)
    {
      for (uint16_t yy = 0; yy < m_Sinogram->N_t; yy++)
      {
        if(m_Sinogram->counts->d[zz][xx][yy] != 0)
        {
#ifdef BRIGHT_FIELD
          Weight->d[zz][xx][yy] = m_Sinogram->counts->d[zz][xx][yy];
#else
          Weight->d[zz][xx][yy] = 1.0 / (m_Sinogram->counts->d[zz][xx][yy] * NuisanceParams.alpha->d[zz]); //The variance for each view ~=averagecounts in that view
#endif
        }
      }
    }
  }
  STOP_TIMER;
  std::cout << std::endl;
  printf("Computing Weights Time ");
  PRINT_TIME;
#endif//Noise model


	for (int16_t i_theta = 0;i_theta < m_Sinogram->N_theta; i_theta++)//slice index
	{
    checksum = 0;
    for (int16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {

        // checksum+=Y_Est->d[i][j][k];
        ErrorSino->d[i_theta][i_r][i_t] = m_Sinogram->counts->d[i_theta][i_r][i_t] - Y_Est->d[i_theta][i_r][i_t] - NuisanceParams.mu->d[i_theta];
        //	if(fabs(ErrorSino->d[i_theta][i_r][i_t]) > 20000)
        //		printf("%d %d %d %lf %lf %lf\n",i_theta,i_r,i_t,ErrorSino->d[i_theta][i_r][i_t],Weight->d[i_theta][i_r][i_t],Y_Est->d[i_theta][i_r][i_t]);

#ifndef NOISE_MODEL
        if(m_Sinogram->counts->d[i_theta][i_r][i_t] != 0)
        {
#ifdef BRIGHT_FIELD
          Weight->d[i_theta][i_r][i_t] = m_Sinogram->counts->d[i_theta][i_r][i_t];
#else
          Weight->d[i_theta][i_r][i_t] = 1.0 / m_Sinogram->counts->d[i_theta][i_r][i_t];
#endif //bright field check
        }
        else Weight->d[i_theta][i_r][i_t] = 0;
#endif//Noise model

#ifdef FORWARD_PROJECT_MODE
        temp=Y_Est->d[i_theta][i_r][i_t]/NuisanceParams.I_0->d[i_theta];
        fwrite(&temp,sizeof(DATA_TYPE),1,Fp6);
#endif
        checksum += Weight->d[i_theta][i_r][i_t];
      }
    }
		printf("Check sum of Diagonal Covariance Matrix= %lf\n", checksum);

  }


#ifdef FORWARD_PROJECT_MODE
	return 0;//exit the program once we finish forward projecting the object
#endif

#ifdef COST_CALCULATE
	DATA_TYPE cost_value;
  std::string filepath(m_Inputs->outputDir);
  filepath = filepath.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::CostFunctionFile);

  CostData::Pointer cost = CostData::New();
  cost->initOutputFile(filepath);

  /*********************Cost Calculation*************************************/
  cost_value = computeCost(ErrorSino, Weight);
  int increase = cost->addCostValue(cost_value);
  if (increase ==1)
  {
	  //Only one cost so dont worry about increase
 //   std::cout << "Cost just increased!" << std::endl;
  }
  cost->writeCostValue(cost_value);
  /**************************************************************************/
#endif //Cost calculation endif


	//Loop through every voxel updating it by solving a cost function
	for(int16_t OuterIter = 0; OuterIter < m_Inputs->NumOuterIter; OuterIter++)
	{
	  indent = "";
	  std::cout << "Outer Iteration " << OuterIter << std::endl;
    for (int16_t Iter = 0; Iter < m_Inputs->NumIter; Iter++)
    {
      indent = "  ";
      std::cout << indent << "Inner Iteration " << Iter << std::endl;
      indent = "    ";
      // This is all done PRIOR to calling what will become a method
      VoxelUpdateType updateType = RegularRandomOrderUpdate;
#ifdef NHICD
      if(0 == Iter % 2) {
        updateType = HomogeniousUpdate;
      }
      else
      {
        updateType = NonHomogeniousUpdate;
      }
#else
		
#endif

      // This could contain multiple Subloops also
 uint8_t status= updateVoxels(Iter, updateType, VisitCount, RandomNumber, TempCol, ErrorSino, Weight, VoxelLineResponse, NuisanceParams, Mask, cost);

 if(status == 0)
   {
     break;//stop inner loop if we have hit the threshold value for x
   }
      // Check to see if we are canceled.
      if(getCancel() == true)
      {
        setErrorCondition(-100);
        return;
      }


    } /* ++++++++++ END Inner Iteration Loop +++++++++++++++ */
#ifdef JOINT_ESTIMATION

    //high=5e100;//this maintains the max and min bracket values for rooting lambda
    DATA_TYPE high = std::numeric_limits<DATA_TYPE>::max();
    DATA_TYPE low = std::numeric_limits<DATA_TYPE>::min();
    //Joint Scale And Offset Estimation

    //forward project
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          Y_Est->d[i_theta][i_r][i_t] = 0.0;
        }
      }
    }

    for (j = 0; j < m_Geometry->N_z; j++)
    {
      for (k = 0; k < m_Geometry->N_x; k++)
      {
        if(TempCol[j][k]->count > 0)
        {
          for (i = 0; i < m_Geometry->N_y; i++) //slice index
          {

            for (uint32_t q = 0; q < TempCol[j][k]->count; q++)
            {
              //calculating the footprint of the voxel in the t-direction

              uint16_t i_theta = int16_t(floor(static_cast<float>(TempCol[j][k]->index[q] / (m_Sinogram->N_r))));
              uint16_t i_r = (TempCol[j][k]->index[q] % (m_Sinogram->N_r));

              uint16_t VoxelLineAccessCounter = 0;
              for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
              {
                Y_Est->d[i_theta][i_r][i_t] += ((TempCol[j][k]->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter++] * m_Geometry->Object->d[j][k][i]));
              }
            }

          }
        }

      }
    }

    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      DATA_TYPE a = 0;
      DATA_TYPE b = 0;
      DATA_TYPE c = 0;
      DATA_TYPE d = 0;
      DATA_TYPE e = 0;
      DATA_TYPE numerator_sum = 0;
      DATA_TYPE denominator_sum = 0;

      //compute the parameters of the quadratic for each view
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {

          numerator_sum += (m_Sinogram->counts->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t]);
          denominator_sum += (Weight->d[i_theta][i_r][i_t]);

          a += (Y_Est->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t]);
          b += (Y_Est->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t] * m_Sinogram->counts->d[i_theta][i_r][i_t]);
          c += (m_Sinogram->counts->d[i_theta][i_r][i_t] * m_Sinogram->counts->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t]);
          d += (Y_Est->d[i_theta][i_r][i_t] * Y_Est->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t]);

        }
      }

      bk_cost->d[i_theta][1] = numerator_sum; //yt*\lambda*1
      bk_cost->d[i_theta][0] = b;//yt*\lambda*(Ax)
      ck_cost->d[i_theta] = c;//yt*\lambda*y
      Qk_cost->d[i_theta][2] = denominator_sum;
      Qk_cost->d[i_theta][1] = a;
      Qk_cost->d[i_theta][0] = d;

      d1->d[i_theta] = numerator_sum / denominator_sum;
      d2->d[i_theta] = a / denominator_sum;

      a = 0;
      b = 0;
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          a += ((Y_Est->d[i_theta][i_r][i_t] - d2->d[i_theta]) * Weight->d[i_theta][i_r][i_t] * Y_Est->d[i_theta][i_r][i_t]);
          b -= ((m_Sinogram->counts->d[i_theta][i_r][i_t] - d1->d[i_theta]) * Weight->d[i_theta][i_r][i_t] * Y_Est->d[i_theta][i_r][i_t]);
        }
      }
      QuadraticParameters->d[i_theta][0] = a;
      QuadraticParameters->d[i_theta][1] = b;

      temp = (QuadraticParameters->d[i_theta][1] * QuadraticParameters->d[i_theta][1]) / (4 * QuadraticParameters->d[i_theta][0]);

      if(temp > 0 && temp < high) high = temp; //high holds the maximum value for the rooting operation. beyond this value discriminants become negative. Basically high = min{b^2/4*a}
      else if(temp < 0 && temp > low) low = temp;

    }
    //compute cost
    /********************************************************************************************/
    sum = 0;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      sum += (Qk_cost->d[i_theta][0] * NuisanceParams.I_0->d[i_theta] * NuisanceParams.I_0->d[i_theta]
          + 2 * Qk_cost->d[i_theta][1] * NuisanceParams.I_0->d[i_theta] * NuisanceParams.mu->d[i_theta]
          + NuisanceParams.mu->d[i_theta] * NuisanceParams.mu->d[i_theta] * Qk_cost->d[i_theta][2]
          - 2 * (bk_cost->d[i_theta][0] * NuisanceParams.I_0->d[i_theta] + NuisanceParams.mu->d[i_theta] * bk_cost->d[i_theta][1]) + ck_cost->d[i_theta]); //evaluating the cost function
    }
    sum /= 2;
    printf("The value of the data match error prior to updating the I and mu =%lf\n", sum);

    /********************************************************************************************/

#ifdef GEOMETRIC_MEAN_CONSTRAINT
    calculateGeometricMeanConstraint();
#else
    DATA_TYPE sum1 = 0;
    DATA_TYPE sum2 = 0;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      sum1 += (1.0 / (Qk_cost->d[i_theta][0] - Qk_cost->d[i_theta][1] * d2->d[i_theta]));
      sum2 += ((bk_cost->d[i_theta][0] - Qk_cost->d[i_theta][1] * d1->d[i_theta]) / (Qk_cost->d[i_theta][0] - Qk_cost->d[i_theta][1] * d2->d[i_theta]));
    }
    DATA_TYPE LagrangeMultiplier = (-m_Sinogram->N_theta * m_Sinogram->TargetGain + sum2) / sum1;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {

      NuisanceParams.I_0->d[i_theta] = (-1 * LagrangeMultiplier - Qk_cost->d[i_theta][1] * d1->d[i_theta] + bk_cost->d[i_theta][0])
      / (Qk_cost->d[i_theta][0] - Qk_cost->d[i_theta][1] * d2->d[i_theta]);

      //Postivity Constraing on the gains

      if (NuisanceParams.I_0->d[i_theta] < 0)
      {
        NuisanceParams.I_0->d[i_theta]=0;
      }

      NuisanceParams.mu->d[i_theta] = d1->d[i_theta] - d2->d[i_theta] * NuisanceParams.I_0->d[i_theta]; //some function of I_0[i_theta]
    }

#endif //Type of constraing Geometric or arithmetic
    //Re normalization

    /*sum=0;
     for(i_theta = 0 ; i_theta < Sinogram->N_theta; i_theta++)
     sum+=(log(NuisanceParams.I_0->d[i_theta]));
     sum/=Sinogram->N_theta;
     temp=exp(sum) - Sinogram->TargetGain;

     printf("%lf\n",temp);*/

    /********************************************************************************************/
    //checking to see if the cost went down
    sum = 0;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      sum += (Qk_cost->d[i_theta][0] * NuisanceParams.I_0->d[i_theta] * NuisanceParams.I_0->d[i_theta]
          + 2 * Qk_cost->d[i_theta][1] * NuisanceParams.I_0->d[i_theta] * NuisanceParams.mu->d[i_theta]
          + NuisanceParams.mu->d[i_theta] * NuisanceParams.mu->d[i_theta] * Qk_cost->d[i_theta][2]
          - 2 * (bk_cost->d[i_theta][0] * NuisanceParams.I_0->d[i_theta] + NuisanceParams.mu->d[i_theta] * bk_cost->d[i_theta][1]) + ck_cost->d[i_theta]); //evaluating the cost function
    }
    sum /= 2;

    printf("The value of the data match error after updating the I and mu =%lf\n", sum);

    /*****************************************************************************************************/

    //Reproject to compute Error Sinogram for ICD
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          ErrorSino->d[i_theta][i_r][i_t] = m_Sinogram->counts->d[i_theta][i_r][i_t] - NuisanceParams.mu->d[i_theta]
          - (NuisanceParams.I_0->d[i_theta] * Y_Est->d[i_theta][i_r][i_t]);
        }
      }
    }

    //Printing the gain and offset after updating
    /*
     printf("Offsets\n");
     for(i_theta = 0 ; i_theta < Sinogram->N_theta; i_theta++)
     {
     printf("%lf\n",NuisanceParams.mu->d[i_theta]);
     }
     printf("Gain\n");
     for(i_theta = 0 ; i_theta < Sinogram->N_theta; i_theta++)
     {
     printf("%lf\n",NuisanceParams.I_0->d[i_theta]);
     }*/

#ifdef COST_CALCULATE
		
		/*********************Cost Calculation*************************************/
		DATA_TYPE cost_value = computeCost(ErrorSino, Weight);
		int increase = cost->addCostValue(cost_value);
		if (increase ==1)
		{
			std::cout << "Cost just increased!" << std::endl;
			break;
		}
		cost->writeCostValue(cost_value);
		/**************************************************************************/
		
#endif
    printf("Lagrange Multiplier = %lf\n", LagrangeMultiplier);

    printf("Offsets\n");
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      printf("%lf\n", NuisanceParams.mu->d[i_theta]);
    }
    printf("Gain\n");
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      printf("%lf\n", NuisanceParams.I_0->d[i_theta]);
    }

#ifdef NOISE_MODEL
    //Updating the Weights
    for(uint16_t i_theta=0;i_theta < m_Sinogram->N_theta;i_theta++)
    {
      sum=0;

      for(uint16_t i_r=0; i_r < m_Sinogram->N_r; i_r++)
      for(uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      sum+=(ErrorSino->d[i_theta][i_r][i_t]*ErrorSino->d[i_theta][i_r][i_t]*Weight->d[i_theta][i_r][i_t]);
      sum/=(m_Sinogram->N_r*m_Sinogram->N_t);

      NuisanceParams.alpha->d[i_theta]=sum;

      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {

          if(NuisanceParams.alpha->d[i_theta] != 0)
          {
            Weight->d[i_theta][i_r][i_t] /= NuisanceParams.alpha->d[i_theta];
          }
        }
      }

    }

#ifdef COST_CALCULATE
		
		/*********************Cost Calculation*************************************/
		cost_value = computeCost(ErrorSino, Weight);
		increase = cost->addCostValue(cost_value);
		if (increase ==1)
		{
			std::cout << "Cost just increased!" << std::endl;
			break;
		}
		cost->writeCostValue(cost_value);
		/**************************************************************************/
		
#endif//cost
#endif//NOISE_MODEL
#endif//Joint estimation endif

  }/* ++++++++++ END Outer Iteration Loop +++++++++++++++ */

  indent = "";
  cost->printCosts(std::cout);


#ifdef FORWARD_PROJECT_MODE
  MAKE_OUTPUT_FILE(Fp6, fileError, m_Inputs->outputDir, ScaleOffsetCorrection::ForwardProjectedObjectFile);
  if (fileError == 1)
  {

  }
#endif


#ifdef DEBUG_CONSTRAINT_OPT
  FILE *Fp8 = NULL;
  MAKE_OUTPUT_FILE(Fp8, fileError, m_Inputs->outputDir, ScaleOffsetCorrection::CostFunctionCoefficientsFile);
  if (fileError >= 0)
  {
    // Write the output File
  }
#endif

/* Write the Gains and Offsets to an output file */
  NuisanceParamWriter::Pointer nuisanceBinWriter = NuisanceParamWriter::New();
  nuisanceBinWriter->setSinogram(m_Sinogram);
  nuisanceBinWriter->setInputs(m_Inputs);
  nuisanceBinWriter->addObserver(this);
  nuisanceBinWriter->setNuisanceParams(&NuisanceParams);
  nuisanceBinWriter->setFileName(ScaleOffsetCorrection::FinalGainParametersFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_I_O);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    updateProgressAndMessage(nuisanceBinWriter->getErrorMessage().c_str(), 100);
  }

  nuisanceBinWriter->setFileName(ScaleOffsetCorrection::FinalOffsetParametersFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_mu);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
   setErrorCondition(-1);
   updateProgressAndMessage(nuisanceBinWriter->getErrorMessage().c_str(), 100);
  }

#ifdef NOISE_MODEL
  nuisanceBinWriter->setFileName(ScaleOffsetCorrection::FinalVariancesFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_alpha);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    updateProgressAndMessage(nuisanceBinWriter->getErrorMessage().c_str(), 100);
  }
#endif//Noise Model


  START_TIMER;
  //calculates Ax and returns a pointer to the memory block
  Final_Sinogram=forwardProject(detectorResponse, H_t);
  STOP_TIMER;
  PRINT_TIME;

  // Write the Sinogram out to a file
  SinogramBinWriter::Pointer sinogramWriter = SinogramBinWriter::New();
  sinogramWriter->setSinogram(m_Sinogram);
  sinogramWriter->setInputs(m_Inputs);
  sinogramWriter->addObserver(this);
  sinogramWriter->setNuisanceParams(&NuisanceParams);
  sinogramWriter->setData(Final_Sinogram);
  sinogramWriter->execute();
  if (sinogramWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    updateProgressAndMessage(sinogramWriter->getErrorMessage().c_str(), 100);
  }

  // Write the Reconstruction out to a file
  RawGeometryWriter::Pointer writer = RawGeometryWriter::New();
  writer->setGeometry(m_Geometry);
  writer->setFilePath(m_Inputs->OutputFile);
  writer->addObserver(this);
  writer->execute();
  if (writer->getErrorCondition() < 0)
  {
    setErrorCondition(writer->getErrorCondition());
    updateProgressAndMessage("Error Writing the Raw Geometry", 100);
  }



  std::string vtkFile(m_Inputs->outputDir);
  vtkFile = vtkFile.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::VtkGeometryFile);

  VTKStructuredPointsFileWriter vtkWriter;
  vtkWriter.setWriteBinaryFiles(true);
  DimsAndRes dimsAndRes;
  dimsAndRes.dim0 = m_Geometry->N_x;
  dimsAndRes.dim1 = m_Geometry->N_y;
  dimsAndRes.dim2 = m_Geometry->N_z;
  dimsAndRes.resx = 1.0f;
  dimsAndRes.resy = 1.0f;
  dimsAndRes.resz = 1.0f;

  std::vector<VtkScalarWriter*> scalarsToWrite;

  VtkScalarWriter* w0 = static_cast<VtkScalarWriter*>(new TomoOutputScalarWriter(m_Geometry));
  w0->setWriteBinaryFiles(true);
  scalarsToWrite.push_back(w0);

  int error = vtkWriter.write<DimsAndRes>(vtkFile, &dimsAndRes, scalarsToWrite);
  if (error < 0)
  {
    std::cout << "Error writing vtk file '" << vtkFile << "'" << std::endl;
  }

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


  updateProgressAndMessage("Reconstruction Complete", 100);
  setErrorCondition(0);
  std::cout << "Total Running Time for Execute: " << (EIMTOMO_getMilliSeconds()-totalTime)/1000 << std::endl;
  return;
}


/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void SOCEngine::minMax(DATA_TYPE *low,DATA_TYPE *high)
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
			//	printf("\n");
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
RealImageType::Pointer SOCEngine::calculateVoxelProfile()
{
	DATA_TYPE angle,MaxValLineIntegral;
	DATA_TYPE temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;
//	DATA_TYPE** VoxProfile = (DATA_TYPE**)multialloc(sizeof(DATA_TYPE),2,m_Sinogram->N_theta,PROFILE_RESOLUTION);
	size_t dims[2] = {m_Sinogram->N_theta,PROFILE_RESOLUTION};
	RealImageType::Pointer VoxProfile = RealImageType::New(dims);
	VoxProfile->setName("VoxelProfile");

	DATA_TYPE checksum=0;
	uint16_t i,j;
	FILE* Fp = NULL;
	int fileError = 0;
  MAKE_OUTPUT_FILE(Fp, fileError, m_Inputs->outputDir, ScaleOffsetCorrection::VoxelProfileFile);
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
			MaxValLineIntegral = m_Inputs->delta_xz/cos(angle);
		}
		else
		{
			MaxValLineIntegral = m_Inputs->delta_xz/cos(M_PI_2-angle);
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
				VoxProfile->d[i][j] = 0;
			else if(t > RightNear)
				VoxProfile->d[i][j] = MaxValLineIntegral*(RightCorner-t)/(RightCorner-RightNear);
			else if(t >= LeftNear)
				VoxProfile->d[i][j] = MaxValLineIntegral;
			else
				VoxProfile->d[i][j] = MaxValLineIntegral*(t-LeftCorner)/(LeftNear-LeftCorner);

			fwrite(&(VoxProfile->d[i][j]),sizeof(DATA_TYPE),1,Fp);
			checksum+=VoxProfile->d[i][j];
		}

	}

	//printf("Pixel Profile Check sum =%lf\n",checksum);
	fclose(Fp);
	return VoxProfile;
}

/*******************************************************************
 Forwards Projects the Object and stores it in a 3-D matrix
 ********************************************************************/
RealVolumeType::Pointer SOCEngine::forwardProject(RealVolumeType::Pointer DetectorResponse, RealVolumeType::Pointer H_t)
{
  updateProgressAndMessage("Executing Forward Projection", 50);

	DATA_TYPE x,z,y;
	DATA_TYPE r,rmin,rmax,t,tmin,tmax;
	DATA_TYPE w1,w2,f1,f2;
	DATA_TYPE center_r,delta_r,center_t,delta_t;
	int16_t index_min,index_max,slice_index_min,slice_index_max,i_r,i_t,i_theta;
	int16_t i,j,k;
	int16_t index_delta_t,index_delta_r;
	RealVolumeType::Pointer Y_Est;
	size_t dims[3] =
  { m_Sinogram->N_theta, m_Sinogram->N_r, m_Sinogram->N_t };
	Y_Est = RealVolumeType::New(dims);
	Y_Est->setName("Y_Est");
//	Y_Est=(RealVolumeType::Pointer )get_3D(m_Sinogram->N_theta,m_Sinogram->N_r,m_Sinogram->N_t,sizeof(DATA_TYPE));

	for (j = 0; j < m_Geometry->N_z; j++)
	{
		for (k=0; k < m_Geometry->N_x; k++)
		{
			x = m_Geometry->x0 + ((DATA_TYPE)k+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
			z = m_Geometry->z0 + ((DATA_TYPE)j+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. z_0 is the left corner

			for (i = 0; i < m_Geometry->N_y ; i++)
			{
				for(i_theta = 0;i_theta < m_Sinogram->N_theta;i_theta++)
				{
					r = x*cosine->d[i_theta] - z*sine->d[i_theta];
					rmin = r - m_Inputs->delta_xz;
					rmax = r + m_Inputs->delta_xz;

					if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax)
						continue;

					index_min = floor(((rmin - m_Sinogram->R0)/m_Sinogram->delta_r));
					index_max = floor((rmax - m_Sinogram->R0)/m_Sinogram->delta_r);

					if(index_max >= m_Sinogram->N_r)
						index_max = m_Sinogram->N_r - 1;

					if(index_min < 0)
						index_min = 0;

					y = m_Geometry->y0 + ((double)i+ 0.5)*m_Inputs->delta_xy;
					t = y;


					tmin = (t - m_Inputs->delta_xy/2) > m_Sinogram->T0 ? t-m_Inputs->delta_xy/2 : m_Sinogram->T0;
					tmax = (t + m_Inputs->delta_xy/2) <= m_Sinogram->TMax? t + m_Inputs->delta_xy/2 : m_Sinogram->TMax;

					slice_index_min = floor((tmin - m_Sinogram->T0)/m_Sinogram->delta_t);
					slice_index_max = floor((tmax - m_Sinogram->T0)/m_Sinogram->delta_t);

					if(slice_index_min < 0)
						slice_index_min = 0;
					if(slice_index_max >= m_Sinogram->N_t)
						slice_index_max = m_Sinogram->N_t-1;

					for (i_r = index_min; i_r <= index_max; i_r++)
					{
						center_r = m_Sinogram->R0 + ((double)i_r + 0.5)*m_Sinogram->delta_r ;
						delta_r = fabs(center_r - r);
						index_delta_r = floor((delta_r/OffsetR));

						if(index_delta_r < DETECTOR_RESPONSE_BINS)
						{
							w1 = delta_r - index_delta_r*OffsetR;
							w2 = (index_delta_r+1)*OffsetR - delta_r;
							f1 = (w2/OffsetR)*DetectorResponse->d[0][i_theta][index_delta_r] + (w1/OffsetR)*DetectorResponse->d[0][i_theta][index_delta_r+1 < DETECTOR_RESPONSE_BINS ? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];
						}
						else
						{
							f1=0;
						}



						for (i_t = slice_index_min; i_t <= slice_index_max; i_t ++)
						{
							center_t = m_Sinogram->T0 + ((double)i_t + 0.5)*m_Sinogram->delta_t;
							delta_t = fabs(center_t - t);
							index_delta_t = floor((delta_t/OffsetT));

							if(index_delta_t < DETECTOR_RESPONSE_BINS)
							{
								w1 = delta_t - index_delta_t*OffsetT;
								w2 = (index_delta_t+1)*OffsetT - delta_t;
								f2 = (w2/OffsetT)*H_t->d[0][i_theta][index_delta_t] + (w1/OffsetT)*H_t->d[0][i_theta][index_delta_t+1 < DETECTOR_RESPONSE_BINS ? index_delta_t+1:DETECTOR_RESPONSE_BINS-1];
							}
							else
							{
								f2 = 0;
							}

							Y_Est->d[i_theta][i_r][i_t]+=(f1*f2*m_Geometry->Object->d[j][k][i]);

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
	cosine = RealArrayType::New(dims);
	cosine->setName("cosine");
	sine = RealArrayType::New(dims);
	sine->setName("sine");

	for(i=0;i<m_Sinogram->N_theta;i++)
	{
		cosine->d[i]=cos(m_Sinogram->angles[i]);
		sine->d[i]=sin(m_Sinogram->angles[i]);
	}
}

void SOCEngine::initializeBeamProfile()
{
	uint16_t i;
	DATA_TYPE sum=0,W;
//	BeamProfile=(DATA_TYPE*)get_spc(BEAM_RESOLUTION,sizeof(DATA_TYPE));
	size_t dims[1] = { BEAM_RESOLUTION };
	BeamProfile = RealArrayType::New(dims);
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
		BeamProfile->d[i]/=m_Sinogram->delta_t;
		// printf("%lf\n",BeamProfile->d[i]);
	}



}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DATA_TYPE SOCEngine::computeCost(RealVolumeType::Pointer ErrorSino,RealVolumeType::Pointer Weight)
{
	DATA_TYPE cost=0,temp=0;
//	DATA_TYPE delta;
	int16_t i,j,k;
//	int16_t p,q,r;

//Data Mismatch Error
	for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (j = 0; j < m_Sinogram->N_r; j++)
    {
      for (k = 0; k < m_Sinogram->N_t; k++)
      {
        cost += (ErrorSino->d[i][j][k] * ErrorSino->d[i][j][k] * Weight->d[i][j][k]);
      }
    }
  }

  cost /= 2;

//  std::cout << "\nCompute Cost: Data mismatch term = " << cost;
//	fflush(stdout);

//Prior Model Error
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
	for (i = 0; i < m_Geometry->N_z; i++)
		for (j = 0; j < m_Geometry->N_x; j++)
			for(k = 0; k < m_Geometry->N_y; k++)
			{

				if(k+1 <  m_Geometry->N_y)
				{
					delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j][k+1];
					temp += FILTER[2][1][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
				}



				if(j+1 < m_Geometry->N_x)
				{
					if(k-1 >= 0)
					{
						delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k-1];
						temp += FILTER[0][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}

					delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k];
					temp += FILTER[1][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));


					if(k+1 < m_Geometry->N_y)
					{
						delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i][j+1][k+1];
						temp += FILTER[2][1][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}

				}

				if(i+1 < m_Geometry->N_z)
				{

					if(j-1 >= 0)
					{
						delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k];
						temp += FILTER[1][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}

					delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k];
					temp += FILTER[1][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));

					if(j+1 < m_Geometry->N_x)
					{
						delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k];
						temp += FILTER[1][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}


					if(j-1 >= 0)
					{
						if(k-1 >= 0)
						{
							delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k-1];
							temp += FILTER[0][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
						}

						if(k+1 < m_Geometry->N_y)
						{
							delta=m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j-1][k+1];
							temp += FILTER[2][2][0]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
						}

					}

					if(k-1 >= 0)
					{
						delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k-1];
						temp += FILTER[0][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}

					if(j+1 < m_Geometry->N_x)
					{
						if(k-1 >= 0)
						{
							delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k-1];
							temp += FILTER[0][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
						}

						if(k+1 < m_Geometry->N_y)
						{
							delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j+1][k+1];
							temp+= FILTER[2][2][2]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
						}
					}

					if(k+1 < m_Geometry->N_y)
					{
						delta = m_Geometry->Object->d[i][j][k]-m_Geometry->Object->d[i+1][j][k+1];
						temp+= FILTER[2][2][1]*(pow(fabs(delta),MRF_P))/(1+pow(fabs(delta/MRF_C), MRF_P-MRF_Q));
					}
				}
			}
			cost+=(temp);
#endif //QGGMRF

	//printf("Cost calculation End..\n");


//Noise Error
#ifdef NOISE_MODEL
	temp = 0;
  for (i = 0; i < m_Sinogram->N_theta; i++)
  {
    if(Weight->d[i][0][0] != 0)
    {
      temp += log(2 * M_PI * (1.0 / Weight->d[i][0][0])); //2*pi*sigma_k^{2}
    }
  }
  temp *= ((m_Sinogram->N_r * m_Sinogram->N_t) / 2);
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
	DATA_TYPE x,z,y;
	DATA_TYPE r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
	DATA_TYPE t; //this is similar to r but along the y direction
	DATA_TYPE tmin,tmax;
	DATA_TYPE rmax,rmin;//stores the start and end points of the pixel profile on the detector
	DATA_TYPE R_Center,TempConst,checksum = 0,delta_r;
//	DATA_TYPE Integral = 0;
	DATA_TYPE T_Center,delta_t;
	DATA_TYPE MaximumSpacePerColumn;//we will use this to allocate space
	DATA_TYPE AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
//	DATA_TYPE ProfileThickness,stepsize;

	//interpolation variables
	DATA_TYPE w1,w2,w3,w4,f1,InterpolatedValue,ContributionAlongT;
//	DATA_TYPE f2;
	int32_t index_min,index_max,slice_index_min,slice_index_max,index_delta_r,index_delta_t;//stores the detector index in which the profile lies
	int32_t BaseIndex,FinalIndex;
//	int32_t ProfileIndex=0;
//	int32_t NumOfDisplacements=32;
	uint32_t count = 0;

  sliceidx = 0;

	AMatrixCol* Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

	x = m_Geometry->x0 + ((DATA_TYPE)col+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	z = m_Geometry->z0 + ((DATA_TYPE)row+0.5)*m_Inputs->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	y = m_Geometry->y0 + ((DATA_TYPE)slice + 0.5)*m_Inputs->delta_xy;

	TempConst=(PROFILE_RESOLUTION)/(2*m_Inputs->delta_xz);

	//alternately over estimate the maximum size require for a single AMatrix column
	AvgNumXElements = ceil(3*m_Inputs->delta_xz/m_Sinogram->delta_r);
	AvgNumYElements = ceil(3*m_Inputs->delta_xy/m_Sinogram->delta_t);
	MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*m_Sinogram->N_theta;

	Temp->values = (DATA_TYPE*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(DATA_TYPE));
	Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));


#ifdef AREA_WEIGHTED
	for(uint32_t i=0;i<m_Sinogram->N_theta;i++)
	{

		r = x*cosine->d[i] - z*sine->d[i];
		t = y;

		rmin = r - m_Inputs->delta_xz;
		rmax = r + m_Inputs->delta_xz;

		tmin = (t - m_Inputs->delta_xy/2) > m_Sinogram->T0 ? t-m_Inputs->delta_xy/2 : m_Sinogram->T0;
		tmax = (t + m_Inputs->delta_xy/2) <= m_Sinogram->TMax ? t + m_Inputs->delta_xy/2 : m_Sinogram->TMax;

		if(rmax < m_Sinogram->R0 || rmin > m_Sinogram->RMax)
			continue;



		index_min = floor(((rmin - m_Sinogram->R0)/m_Sinogram->delta_r));
		index_max = floor((rmax - m_Sinogram->R0)/m_Sinogram->delta_r);


		if(index_max >= m_Sinogram->N_r)
			index_max = m_Sinogram->N_r - 1;

		if(index_min < 0)
			index_min = 0;

		slice_index_min = floor((tmin - m_Sinogram->T0)/m_Sinogram->delta_t);
		slice_index_max = floor((tmax - m_Sinogram->T0)/m_Sinogram->delta_t);

		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= m_Sinogram->N_t)
			slice_index_max = m_Sinogram->N_t -1;

		BaseIndex = i*m_Sinogram->N_r;//*Sinogram->N_t;

		for(j = index_min;j <= index_max; j++)//Check
		{

			//Accounting for Beam width
			R_Center = (m_Sinogram->R0 + (((DATA_TYPE)j) + 0.5) *(m_Sinogram->delta_r));//the 0.5 is to get to the center of the detector

			//Find the difference between the center of detector and center of projection and compute the Index to look up into
			delta_r = fabs(r - R_Center);
			index_delta_r = floor((delta_r/OffsetR));


			if (index_delta_r >= 0 && index_delta_r < DETECTOR_RESPONSE_BINS)
			{

		//		for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
		//		{
					T_Center = (m_Sinogram->T0 + (((DATA_TYPE)sliceidx) + 0.5) *(m_Sinogram->delta_t));
					delta_t = fabs(t - T_Center);
					index_delta_t = 0;//floor(delta_t/OffsetT);



					if (index_delta_t >= 0 && index_delta_t < DETECTOR_RESPONSE_BINS)
					{

						//Using index_delta_t,index_delta_t+1,index_delta_r and index_delta_r+1 do bilinear interpolation
						w1 = delta_r - index_delta_r*OffsetR;
						w2 = (index_delta_r+1)*OffsetR - delta_r;

						w3 = delta_t - index_delta_t*OffsetT;
						w4 = (index_delta_r+1)*OffsetT - delta_t;


						f1 = (w2/OffsetR)*DetectorResponse->d[index_delta_t][i][index_delta_r] + (w1/OffsetR)*DetectorResponse->d[index_delta_t][i][index_delta_r+1 < DETECTOR_RESPONSE_BINS ? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];
						//	f2 = (w2/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS ?index_delta_t+1 : DETECTOR_RESPONSE_BINS-1][i][index_delta_r] + (w1/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS? index_delta_t+1:DETECTOR_RESPONSE_BINS][i][index_delta_r+1 < DETECTOR_RESPONSE_BINS? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];

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


	Ai->values=(DATA_TYPE*)get_spc(count,sizeof(DATA_TYPE));
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
DATA_TYPE SOCEngine::absMaxArray(std::vector<DATA_TYPE> &Array)
{
  uint16_t i;
  DATA_TYPE max;
  max = fabs(Array[0]);
  for(i =1; i < Array.size();i++)
    if(fabs(Array[i]) > max)
      max=fabs(Array[i]);
  return max;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::calculateGeometricMeanConstraint(ScaleOffsetParams* NuisanceParams)
{

  DATA_TYPE low,high,dist;
  DATA_TYPE perturbation=1e-30;//perturbs the rooting range
  uint16_t rooting_attempt_counter;
  int32_t errorcode=-1;
  DATA_TYPE LagrangeMultiplier;
  DATA_TYPE LambdaRootingAccuracy=1e-10;//accuracy for rooting Lambda
  uint16_t MaxNumRootingAttempts = 1000;//for lambda this corresponds to a distance of 2^1000
  DATA_TYPE* root;
  DATA_TYPE a,b;
  DATA_TYPE temp_mu;

//TODO: Fix these initializations
  high = 0;
  low = 0;
  dist = 0;
  //Root the expression using the derived quadratic parameters. Need to choose min and max values
     printf("Rooting the equation to solve for the optimal Lagrange multiplier\n");

     if(high != std::numeric_limits<DATA_TYPE>::max())
     {
       high-=perturbation; //Since the high value is set to make all discriminants exactly >=0 there are some issues when it is very close due to round off issues. So we get sqrt(-6e-20) for example. So subtract an arbitrary value like 0.5
       //we need to find a window within which we need to root the expression . the upper bound is clear but lower bound we need to look for one
       //low=high;
       dist=-1;
     }
     else if (low != std::numeric_limits<DATA_TYPE>::min())
     {
       low +=perturbation;
       //high=low;
       dist=1;
     }

     /*  fa=CE_ConstraintEquation(high);
      signa=fa>0;
      fb=CE_ConstraintEquation(low);
      signb = fb>0;
      dist=1;
      while (signb != signa )
      {
      low = low - dist;
      fb=CE_ConstraintEquation(low);
      signb = fb>0;
      dist*=2;
      }*/
     rooting_attempt_counter=0;
     errorcode=-1;
     CE_ConstraintEquation ce(NumOfViews, QuadraticParameters, d1, d2, Qk_cost, bk_cost, ck_cost, LogGain);

     while(errorcode != 0 && low <= high && rooting_attempt_counter < MaxNumRootingAttempts) //0 implies the signof the function at low and high is the same
     {

       LagrangeMultiplier=solve<CE_ConstraintEquation>(&ce, low, high, LambdaRootingAccuracy, &errorcode);
       low=low+dist;
       dist*=2;
       rooting_attempt_counter++;
     }

     //Something went wrong and the algorithm was unable to bracket the root within the given interval
     if(rooting_attempt_counter == MaxNumRootingAttempts && errorcode != 0)
     {
       printf("The rooting for lambda was unsuccesful\n");
       printf("Low = %lf High = %lf\n",low,high);
       printf("Quadratic Parameters\n");
       for(int i_theta =0; i_theta < m_Sinogram->N_theta;i_theta++)
       {
         printf("%lf %lf %lf\n",QuadraticParameters->d[i_theta][0],QuadraticParameters->d[i_theta][1],QuadraticParameters->d[i_theta][2]);
       }
 #ifdef DEBUG_CONSTRAINT_OPT

       //for(i_theta =0; i_theta < Sinogram->N_theta;i_theta++)
       //{
       fwrite(&Qk_cost->d[0][0], sizeof(DATA_TYPE), m_Sinogram->N_theta*3, Fp8);
       //}
       //for(i_theta =0; i_theta < Sinogram->N_theta;i_theta++)
       //{
       fwrite(&bk_cost->d[0][0], sizeof(DATA_TYPE), m_Sinogram->N_theta*2, Fp8);
       //}
       fwrite(&ck_cost->d[0], sizeof(DATA_TYPE), m_Sinogram->N_theta, Fp8);
 #endif

       return;
     }

     CE_ConstraintEquation conEqn(NumOfViews, QuadraticParameters, d1, d2, Qk_cost, bk_cost, ck_cost, LogGain);

     //Based on the optimal lambda compute the optimal mu and I0 values
     for (int i_theta =0; i_theta < m_Sinogram->N_theta; i_theta++)
     {

       root = conEqn.CE_RootsOfQuadraticFunction(QuadraticParameters->d[i_theta][0],QuadraticParameters->d[i_theta][1],LagrangeMultiplier); //returns the 2 roots of the quadratic parameterized by a,b,c
       a=root[0];
       b=root[0];
       if(root[0] >= 0 && root[1] >= 0)
       {
         temp_mu = d1->d[i_theta] - root[0]*d2->d[i_theta]; //for a given lambda we can calculate I0(\lambda) and hence mu(lambda)
         a = (Qk_cost->d[i_theta][0]*root[0]*root[0] + 2*Qk_cost->d[i_theta][1]*root[0]*temp_mu + temp_mu*temp_mu*Qk_cost->d[i_theta][2] - 2*(bk_cost->d[i_theta][0]*root[0] + temp_mu*bk_cost->d[i_theta][1]) + ck_cost->d[i_theta]);//evaluating the cost function

         temp_mu = d1->d[i_theta] - root[1]*d2->d[i_theta];//for a given lambda we can calculate I0(\lambda) and hence mu(lambda)
         b = (Qk_cost->d[i_theta][0]*root[1]*root[1] + 2*Qk_cost->d[i_theta][1]*root[1]*temp_mu + temp_mu*temp_mu*Qk_cost->d[i_theta][2] - 2*(bk_cost->d[i_theta][0]*root[1] + temp_mu*bk_cost->d[i_theta][1]) + ck_cost->d[i_theta]);//evaluating the cost function
       }

       if(a == Minimum(a, b))
       NuisanceParams->I_0->d[i_theta] = root[0];
       else
       {
         NuisanceParams->I_0->d[i_theta] = root[1];
       }

       free(root); //freeing the memory holding the roots of the quadratic equation

       NuisanceParams->mu->d[i_theta] = d1->d[i_theta] - d2->d[i_theta]*NuisanceParams->I_0->d[i_theta];//some function of I_0[i_theta]
     }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::ComputeVSC()
{
  DATA_TYPE filter_op = 0;
  int err = 0;
  FILE *Fp = NULL;
  MAKE_OUTPUT_FILE(Fp, err, m_Inputs->outputDir, ScaleOffsetCorrection::MagnitudeMapFile);
  if(err < 0)
  {

  }
  fwrite(&(MagUpdateMap->d[0][0]), m_Geometry->N_x * m_Geometry->N_z, sizeof(DATA_TYPE), Fp);
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
            filter_op += HAMMING_WINDOW[p + 2][q + 2] * MagUpdateMap->d[i + p][j + q];
          }
        }
      }
      FiltMagUpdateMap->d[i][j] = filter_op;
    }
  }


  MAKE_OUTPUT_FILE(Fp, err, m_Inputs->outputDir, ScaleOffsetCorrection::FilteredMagMapFile);
  if(err < 0)
  {

  }
  fwrite(&(FiltMagUpdateMap->d[0][0]), m_Geometry->N_x * m_Geometry->N_z, sizeof(DATA_TYPE), Fp);
  fclose(Fp);
}


//Sort the entries of FiltMagUpdateMap and set the threshold to be ? percentile
DATA_TYPE SOCEngine::SetNonHomThreshold()
{
	size_t dims[2]={m_Geometry->N_z*m_Geometry->N_x,0};
	RealArrayType::Pointer TempMagMap=RealArrayType::New(dims);
	//TempMagMap->setName("TempMagMap");

	uint32_t ArrLength=m_Geometry->N_z*m_Geometry->N_x;
	DATA_TYPE threshold;

	//Copy into a linear list for easier partial sorting
	for (uint16_t i=0; i < m_Geometry->N_z; i++)
	    for (uint16_t j=0; j < m_Geometry->N_x; j++)
		{
			TempMagMap->d[i*m_Geometry->N_x+j]=FiltMagUpdateMap->d[i][j];
		}


	//Insertion sort

	int16_t j;
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




uint8_t SOCEngine::updateVoxels(int16_t Iter,
                             VoxelUpdateType updateType,
                             UInt8ImageType::Pointer VisitCount,
                             RNGVars* RandomNumber,
                             AMatrixCol*** TempCol,
                             RealVolumeType::Pointer ErrorSino,
                             RealVolumeType::Pointer Weight,
                             AMatrixCol* VoxelLineResponse,
                             ScaleOffsetParams &NuisanceParams,
                             UInt8ImageType::Pointer Mask,
                             CostData::Pointer cost)
{
  uint8_t exit_status=1;//Indicates normal exit ; else indicates to stop inner iterations
  uint16_t subIterations = 1;
  std::string indent("    ");
  uint8_t err = 0;
  DATA_TYPE UpdatedVoxelValue;
  DATA_TYPE accuracy = 1e-7; //This is the rooting accuracy for x
  int32_t errorcode = -1;
  int16_t Idx;

  //FIXME: Where are these Initialized? Or what values should they be initialized to?
	DATA_TYPE low = 0.0, high = 0.0;

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

  DATA_TYPE NH_Threshold = 0.0;

  for (uint16_t NH_Iter = 0; NH_Iter < subIterations; ++NH_Iter)
  {
    if(updateType == NonHomogeniousUpdate)
    {
      //Compute VSC and create a map of pixels that are above the threshold value
      ComputeVSC();
      NH_Threshold = SetNonHomThreshold();
	  std::cout<<"NHICD THreshold"<<NH_Threshold<<std::endl;
      //Use  FiltMagUpdateMap  to find MagnitudeUpdateMask
      //std::cout << "Completed Calculation of filtered magnitude" << std::endl;
      //Calculate the threshold for the top ? % of voxel updates
    }

    //printf("Iter %d\n",Iter);
#ifdef ROI
    //variables used to stop the process
    DATA_TYPE AverageUpdate = 0;
    DATA_TYPE AverageMagnitudeOfRecon = 0;
#endif

#ifdef RANDOM_ORDER_UPDATES
    int32_t ArraySize = m_Geometry->N_x * m_Geometry->N_z;
    size_t dims[3] =
    { ArraySize, 0, 0 };
    Int32ArrayType::Pointer Counter = Int32ArrayType::New(dims);

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
          if(MagUpdateMap->d[j][k] > NH_Threshold)
          {
            MagUpdateMask->d[j][j] = 1;
            MagUpdateMap->d[j][k] = 0;
	    NumVoxelsToUpdate++;
          }
          else
          {
            MagUpdateMask->d[j][j] = 0;
          }
        }
        else if(updateType == HomogeniousUpdate)
        {
          MagUpdateMap->d[j][k] = 0;
        }
	else if(updateType == RegularRandomOrderUpdate)
	  {
	    MagUpdateMap->d[j][k]=0;
	    NumVoxelsToUpdate++;
	  }

        VisitCount->d[j][k] = 0;
      }
    }
     std::cout<<"Number of voxels to update"<<NumVoxelsToUpdate<<std::endl;
#endif

    START_TIMER;
    for (int32_t j = 0; j < m_Geometry->N_z; j++) //Row index
    {
      for (int32_t k = 0; k < m_Geometry->N_x; k++) //Column index
      {

#ifdef RANDOM_ORDER_UPDATES
        //RandomNumber=init_genrand(Iter);
        int32_t Index = (genrand_int31(RandomNumber)) % ArraySize;
        int32_t k_new = Counter->d[Index] % m_Geometry->N_x;
        int32_t j_new = Counter->d[Index] / m_Geometry->N_x;
        //memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index-1));
        //TODO: Instead just swap the value in Index with the one in ArraySize
        Counter->d[Index] = Counter->d[ArraySize - 1];
        VisitCount->d[j_new][k_new] = 1;
        ArraySize--;
#else
        j_new=j;
        k_new=k;
#endif //Random order updates
        AMatrixCol* TempMemBlock = TempCol[j_new][k_new]; //Remove this

        int shouldInitNeighborhood = 0;
        
        if(updateType == NonHomogeniousUpdate && MagUpdateMask->d[j_new][k_new] == 1 && TempMemBlock->count > 0)
        {
	      ++shouldInitNeighborhood;
        }
	if(updateType == HomogeniousUpdate  && TempMemBlock->count > 0)
	{
	      ++shouldInitNeighborhood;
	}
	if(updateType == RegularRandomOrderUpdate && TempMemBlock->count > 0)
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
            for(q = -1; q <= 1; q++)
            for(r = -1; r <= 1;r++)
            {
              tempindex_x = mod(r+k_new,m_Geometry->N_x);
              tempindex_y =mod(p+i,m_Geometry->N_y);
              tempindex_z = mod(q+j_new,m_Geometry->N_z);
              NEIGHBORHOOD[p+1][q+1][r+1] = m_Geometry->Object->d[tempindex_z][tempindex_x][tempindex_y];
              BOUNDARYFLAG[p+1][q+1][r+1]=1;
            }
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
                        NEIGHBORHOOD[p + 1][q + 1][r + 1] = m_Geometry->Object->d[q + j_new][r + k_new][p + i];
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
              printf("Geom %lf\n", m_Geometry->Object->d[i][31][31]);
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
            V = m_Geometry->Object->d[j_new][k_new][i]; //Store the present value of the voxel
            THETA1 = 0.0;
            THETA2 = 0.0;

            //TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
            for (uint32_t q = 0; q < TempMemBlock->count; q++)
            {
              uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (m_Sinogram->N_r)));
              uint16_t i_r = (TempMemBlock->index[q] % (m_Sinogram->N_r));
              uint16_t VoxelLineAccessCounter = 0;
              for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
              {
                THETA2 += ((NuisanceParams.I_0->d[i_theta] * NuisanceParams.I_0->d[i_theta])
                    * (VoxelLineResponse[i].values[VoxelLineAccessCounter] * VoxelLineResponse[i].values[VoxelLineAccessCounter]) * (TempMemBlock->values[q])
                    * (TempMemBlock->values[q]) * Weight->d[i_theta][i_r][i_t]);
                THETA1 += NuisanceParams.I_0->d[i_theta] * ErrorSino->d[i_theta][i_r][i_t] * (TempMemBlock->values[q])
                    * (VoxelLineResponse[i].values[VoxelLineAccessCounter]) * Weight->d[i_theta][i_r][i_t];
                VoxelLineAccessCounter++;
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

            UpdatedVoxelValue = (DATA_TYPE)solve<DerivOfCostFunc>(&docf, (double)low, (double)high, (double)accuracy, &errorcode);

#else

            errorcode=0;
#ifdef QGGMRF
            UpdatedVoxelValue = CE_FunctionalSubstitution(low,high);

#else
            SurrogateUpdate=surrogateFunctionBasedMin();
            UpdatedVoxelValue=SurrogateUpdate;
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
                printf("Error \n");
                printf("%d %d\n", j_new, k_new);
              }
            }

            //TODO Print appropriate error messages for other values of error code
            m_Geometry->Object->d[j_new][k_new][i] = UpdatedVoxelValue;
	    //#ifdef NHICD
            MagUpdateMap->d[j_new][k_new] += fabs(m_Geometry->Object->d[j_new][k_new][i] - V);
	    //#endif

#ifdef ROI
            if(Mask->d[j_new][k_new] == 1)
            {
              AverageUpdate += fabs(m_Geometry->Object->d[j_new][k_new][i] - V);
              AverageMagnitudeOfRecon += fabs(V); //computing the percentage update =(Change in mag/Initial magnitude)
            }
#endif

            //Update the ErrorSinogram
            for (uint32_t q = 0; q < TempMemBlock->count; q++)
            {
              uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (m_Sinogram->N_r)));
              uint16_t i_r = (TempMemBlock->index[q] % (m_Sinogram->N_r));
              uint16_t VoxelLineAccessCounter = 0;
              for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
              //for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
              {
                ErrorSino->d[i_theta][i_r][i_t] -= (NuisanceParams.I_0->d[i_theta]
                    * (TempMemBlock->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (m_Geometry->Object->d[j_new][k_new][i] - V)));
                VoxelLineAccessCounter++;
              }
            }
            Idx++;
          }
        }
        else
        {
          continue;
        }
      }
    }
    STOP_TIMER;
    PRINT_TIME;

#ifdef RANDOM_ORDER_UPDATES
    for (int j = 0; j < m_Geometry->N_z; j++)
    { //Row index
      for (int k = 0; k < m_Geometry->N_x; k++)
      { //Column index
        if(VisitCount->d[j][k] == 0)
        {
          printf("Pixel (%d %d) not visited\n", j, k);
        }
      }
    }
#endif

#ifdef COST_CALCULATE

    /*********************Cost Calculation*************************************/
    DATA_TYPE cost_value = computeCost(ErrorSino, Weight);
    int increase = cost->addCostValue(cost_value);
    if (increase ==1)
    {
      std::cout << "Cost just increased!" << std::endl;
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
      if((AverageUpdate / AverageMagnitudeOfRecon) < m_Inputs->StopThreshold)
      {
        printf("This is the terminating point %d\n", Iter);
        m_Inputs->StopThreshold *= THRESHOLD_REDUCTION_FACTOR; //Reducing the thresold for subsequent iterations
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
      NumOfBytesWritten=fwrite(&(m_Geometry->Object->d[0][0][0]), sizeof(DATA_TYPE),m_Geometry->N_x*m_Geometry->N_y*m_Geometry->N_z, Fp3);
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
  return exit_status;

}



