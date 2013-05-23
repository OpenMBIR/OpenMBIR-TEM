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

#include "HAADFForwardModel.h"

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
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Common/EIMTime.h"

#include "MBIRLib/GenericFilters/MRCSinogramInitializer.h"

#include "MBIRLib/HAADF/ForwardProject.h"
#include "MBIRLib/HAADF/UpdateYSlice.h"
#include "MBIRLib/HAADF/Filters/NuisanceParamWriter.h"
#include "MBIRLib/HAADF/Filters/NuisanceParamReader.h"
#include "MBIRLib/HAADF/Filters/GainsOffsetsReader.h"
#include "MBIRLib/HAADF/Filters/ComputeInitialOffsets.h"
#include "MBIRLib/HAADF/Filters/SinogramBinWriter.h"

#define USE_TBB_TASK_GROUP 1
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#include <tbb/task.h>
#endif



#define MAKE_OUTPUT_FILE(Fp, outdir, filename)\
    {\
    std::string filepath(outdir);\
    filepath = filepath.append(MXADir::getSeparator()).append(filename);\
    errno = 0;\
    Fp = fopen(filepath.c_str(),"wb");\
    if (Fp == NULL || errno > 0) { std::cout << "Error " << errno << " Opening Output file " << filepath << std::endl;}\
    }

#define START_TIMER uint64_t startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER uint64_t stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME(msg)\
    std::cout << indent << msg << ": " << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFForwardModel::HAADFForwardModel() :
    m_Verbose(false), m_VeryVerbose(false), m_ErrorCondition(0), m_Cancel(false), m_UseDefaultOffset(false)
{


  //Hamming Window here
  k_HammingWindow[0][0] = 0.0013;
  k_HammingWindow[0][1] = 0.0086;
  k_HammingWindow[0][2] = 0.0159;
  k_HammingWindow[0][3] = 0.0086;
  k_HammingWindow[0][4] = 0.0013;
  k_HammingWindow[1][0] = 0.0086;
  k_HammingWindow[1][1] = 0.0581;
  k_HammingWindow[1][2] = 0.1076;
  k_HammingWindow[1][3] = 0.0581;
  k_HammingWindow[1][4] = 0.0086;
  k_HammingWindow[2][0] = 0.0159;
  k_HammingWindow[2][1] = 0.1076;
  k_HammingWindow[2][2] = 0.1993;
  k_HammingWindow[2][3] = 0.1076;
  k_HammingWindow[2][4] = 0.0159;
  k_HammingWindow[3][0] = 0.0013;
  k_HammingWindow[3][1] = 0.0086;
  k_HammingWindow[3][2] = 0.0159;
  k_HammingWindow[3][3] = 0.0086;
  k_HammingWindow[3][4] = 0.0013;
  k_HammingWindow[4][0] = 0.0086;
  k_HammingWindow[4][1] = 0.0581;
  k_HammingWindow[4][2] = 0.1076;
  k_HammingWindow[4][3] = 0.0581;
  k_HammingWindow[4][4] = 0.0086;

  m_TargetGain = 0.0;
  m_InitialGain = RealArrayType::NullPointer();
  m_InitialOffset = RealArrayType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADFForwardModel::~HAADFForwardModel()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::forwardProject(SinogramPtr sinogram,
                                      GeometryPtr geometry,
                                      std::vector<HAADFAMatrixCol::Pointer> &tempCol,
                                      std::vector<HAADFAMatrixCol::Pointer> &voxelLineResponse,
                                      RealVolumeType::Pointer yEstimate,
                                      RealVolumeType::Pointer errorSinogram)
{
  std::string indent("  ");

#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
  tbb::task_scheduler_init init;
  //   int m_NumThreads = init.default_num_threads();
#else
  //   int m_NumThreads = 1;
#endif

  notify("Starting Forward Projection", 10, Observable::UpdateProgressValueAndMessage);
  START_TIMER;
  // This next section looks crazy with all the #if's but this makes sure we are
  // running the exact same code whether in parallel or serial.

#if OpenMBIR_USE_PARALLEL_ALGORITHMS
  tbb::task_group* g = new tbb::task_group;
  if(getVerbose())
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
  for (uint16_t t = 0; t < geometry->N_z; t++)
  {
#if OpenMBIR_USE_PARALLEL_ALGORITHMS
    g->run(ForwardProject(sinogram.get(), geometry.get(), tempCol, voxelLineResponse, yEstimate, this, t, this));
#else
   // ForwardProject fp(sinogram.get(), geometry.get(), tempCol, voxelLineResponse, yEstimate, NuisanceParams.get(), t, this);
   ForwardProject fp(sinogram.get(), geometry.get(), tempCol, voxelLineResponse, yEstimate,this, t, this);
	  //fp.setObservers(getObservers());
    fp();
#endif
  }
#if OpenMBIR_USE_PARALLEL_ALGORITHMS
  g->wait(); // Wait for all the threads to complete before moving on.
  delete g;
#endif

  STOP_TIMER;
  PRINT_TIME("Forward Project Time");

  //Calculate Error Sinogram - Can this be combined with previous loop?
  //Also compute weights of the diagonal covariance matrix
  //Calculate Error Sinogram - Can this be combined with previous loop?
  //Also compute weights of the diagonal covariance matrix
  calculateMeasurementWeight(sinogram, errorSinogram, yEstimate);
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::setQGGMRFValues(QGGMRF::QGGMRF_Values* qggmrf_values)
{
  m_QGGMRF_Values = qggmrf_values;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::printNuisanceParameters(SinogramPtr sinogram)
{
  if(getVeryVerbose())
  {
    // Print out the Initial Gains, Offsets, Variances
    std::cout << "---------------- Initial Gains, Offsets, Variances -------------------" << std::endl;
    std::cout << "Tilt\tGain\tOffset";

    if(NULL != m_InitialVariance.get())
    {
      std::cout << "\tVariance";
    }
    std::cout << std::endl;

    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      std::cout << i_theta << "\t" << m_InitialGain->d[i_theta] << "\t" << m_InitialOffset->d[i_theta];
      if(NULL != m_InitialVariance.get())
      {
        std::cout << "\t" << m_InitialVariance->d[i_theta];
      }
      std::cout << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::resetNuisanceParameters()
{
  m_I_0 = RealArrayType::NullPointer();
  m_Mu = RealArrayType::NullPointer();
  m_Alpha = RealArrayType::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::allocateNuisanceParameters(SinogramPtr sinogram)
{
  //Gain, Offset and Variance Parameter Structures
  size_t dims[3];
  dims[1] = sinogram->N_t;
  dims[0] = sinogram->N_theta;
  m_I_0 = RealArrayType::New(dims, "HAADF::NuisanceParams->I_0");
  m_Mu = RealArrayType::New(dims, "HAADF::NuisanceParams->mu");
  if(m_AdvParams->NOISE_ESTIMATION)
  {
    //alpha is the noise variance adjustment factor
    m_Alpha = RealArrayType::New(dims, "HAADF::NuisanceParams->alpha");
  }
  else
  {
    m_Alpha = RealArrayType::NullPointer();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::gainAndOffsetInitialization(uint16_t nTheta)
{
  Real_t sum = 0;
  Real_t temp = 0;
  for (uint16_t k = 0; k < nTheta; k++)
  {
    // Gains
    m_I_0->d[k] = m_InitialGain->d[k];
    // Offsets
    m_Mu->d[k] = m_InitialOffset->d[k];

    sum += m_I_0->d[k];

  }
  sum /= nTheta;

  if(getVerbose())
  {
    printf("The Arithmetic mean of the constraint is %lf\n", sum);
  }
  if(sum - m_TargetGain > 1e-5)
  {
    if(getVerbose())
    {
      printf("Arithmetic Mean Constraint not met..renormalizing\n");
    }
    temp = m_TargetGain / sum;
    for (uint16_t k = 0; k < nTheta; k++)
    {
      m_I_0->d[k] = m_InitialGain->d[k] * temp;
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::weightInitialization(size_t dims[3])
{
  m_Weight = RealVolumeType::New(dims, "Weight");
#ifdef BF_RECON
	//This variable selects which entries to retain in the sinogram
  m_Selector = UInt8VolumeType::New(dims, "Selector"); 
#endif //BF RECON
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::costInitialization(SinogramPtr sinogram)
{
  size_t dims[3];

  dims[0] = sinogram->N_theta;
  dims[1] = 3;
  dims[2] = 0;

  m_QuadraticParameters = RealImageType::New(dims, "QuadraticParameters");

  m_QkCost = RealImageType::New(dims, "Qk_cost");
  dims[1] = 2;
  m_BkCost = RealImageType::New(dims, "bk_cost");

  dims[0] = sinogram->N_theta;
  m_CkCost = RealArrayType::New(dims, "ck_cost");

  dims[0] = sinogram->N_theta;
  m_D1 = RealArrayType::New(dims, "d1");
  m_D2 = RealArrayType::New(dims, "d2");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::initializeBrightFieldData(SinogramPtr sinogram)
{
  std::stringstream ss;
  if(m_BFTomoInputs.get() != NULL && m_BFSinogram.get() != NULL && m_BFTomoInputs->sinoFile.empty() == false)
  {
    ss << "Initializing BF data";
    notify(ss.str(), 0, Observable::UpdateProgressMessage);

    TomoFilter::Pointer dataReader = TomoFilter::NullPointer();
    std::string extension = MXAFileInfo::extension(m_BFTomoInputs->sinoFile);
    if(extension.compare("mrc") == 0 || extension.compare("ali") == 0)
    {
      dataReader = MRCSinogramInitializer::NewTomoFilter();
    }
    else
    {
      setErrorCondition(-1);
      notify("A supported file reader for the Bright Field file was not found.", 100, Observable::UpdateProgressValueAndMessage);
      return -1;
    }
    dataReader->setTomoInputs(m_BFTomoInputs);
    dataReader->setSinogram(m_BFSinogram);
    dataReader->setAdvParams(m_AdvParams);
    dataReader->setObservers(getObservers());
    dataReader->execute();
    if(dataReader->getErrorCondition() < 0)
    {
      notify("Error reading Input Sinogram Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(dataReader->getErrorCondition());
      return -1;
    }

    //Normalize the HAADF image
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    { //slice index
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
          //1000 is for Marc De Graef data which needed to multiplied
          //  Real_t ttmp = (m_BFSinogram->counts->getValue(i_theta, i_r, i_t) * 1000);
          //  sinogram->counts->divideByValue(ttmp, i_theta, i_r, i_t);
          //100 is for Marc De Graef data which needed to multiplied
          m_BFSinogram->counts->multiplyByValue(100, i_theta, i_r, i_t);
        }
      }
    }

    m_BF_Flag = true;
    notify("BF initialization complete", 0, Observable::UpdateProgressMessage);
  }
  else
  {
    m_BF_Flag = false;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Calculate Error Sinogram
// Also compute weights of the diagonal covariance matrix
// -----------------------------------------------------------------------------
void HAADFForwardModel::calculateMeasurementWeight(SinogramPtr sinogram, RealVolumeType::Pointer errorSinogram, RealVolumeType::Pointer yEstimate)
{
  std::string indent("  ");
  Real_t checksum = 0;
  START_TIMER;
  for (int16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++) //slice index
  {
    if(m_AdvParams->NOISE_ESTIMATION)
    {
      m_Alpha->d[i_theta] = m_InitialVariance->d[i_theta]; //Initialize the refinement parameters from any previous run
    } //Noise model

    checksum = 0;
    for (int16_t i_r = 0; i_r < sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
      {
        size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
        size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
        size_t yest_idx = yEstimate->calcIndex(i_theta, i_r, i_t);
        size_t error_idx = errorSinogram->calcIndex(i_theta, i_r, i_t);
        if(getBF_Flag() == false)
        {
          errorSinogram->d[error_idx] = sinogram->counts->d[counts_idx] - yEstimate->d[yest_idx] - m_Mu->d[i_theta];
        }
        else
        {
          size_t bfcounts_idx = m_BFSinogram->counts->calcIndex(i_theta, i_r, i_t);

          errorSinogram->d[error_idx] = sinogram->counts->d[counts_idx] - m_BFSinogram->counts->d[bfcounts_idx] * yEstimate->d[yest_idx] - m_Mu->d[i_theta];
        }

#ifndef IDENTITY_NOISE_MODEL
        if(sinogram->counts->d[counts_idx] != 0)
        {
          m_Weight->d[weight_idx] = 1.0 / sinogram->counts->d[counts_idx];
        }
        else
        {
          m_Weight->d[weight_idx] = 1.0; //Set the weight to some small number
          //TODO: Make this something resonable
        }
		
		  //If its a bright field recon just over ride the weights
#ifdef BF_RECON
		  m_Weight->d[weight_idx] = BF_MAX/exp(sinogram->counts->d[counts_idx]);  
#endif //BF_RECON
		  
#else
        m_Weight->d[weight_idx] = 1.0;
#endif //IDENTITY_NOISE_MODEL endif
		  
#ifdef BF_RECON
		m_Selector->d[weight_idx]=1;//By default all enties are chosen
#endif //BF_RECON
		  
#ifdef FORWARD_PROJECT_MODE
        temp=yEstimate->d[i_theta][i_r][i_t]/m_I_0->d[i_theta];
        fwrite(&temp,sizeof(Real_t),1,Fp6);
#endif
		  
		  
#ifdef DEBUG
        if(m_Weight->d[weight_idx] < 0)
        {
          //  std::cout << sinogram->counts->d[counts_idx] << "    " << m_Alpha->d[i_theta] << std::endl;
        }
#endif//Debug
        if(m_AdvParams->NOISE_ESTIMATION)
        {
          m_Weight->d[weight_idx] /= m_Alpha->d[i_theta];
        } // NOISE_MODEL

        checksum += m_Weight->d[weight_idx];
      }
    }
    if(getVerbose())
    {
      printf("Check sum of Diagonal Covariance Matrix= %lf\n", checksum);
    }
  }
  STOP_TIMER;
  PRINT_TIME("Computing Weights");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::jointEstimation(SinogramPtr sinogram, RealVolumeType::Pointer errorSinogram, RealVolumeType::Pointer yEstimate, CostData::Pointer cost)
{
  std::stringstream ss;
  std::string indent("  ");
#ifndef BF_RECON
  {
    Real_t AverageI_kUpdate = 0; //absolute sum of the gain updates
    Real_t AverageMagI_k = 0; //absolute sum of the initial gains

    Real_t AverageDelta_kUpdate = 0; //absolute sum of the offsets
    Real_t AverageMagDelta_k = 0; //abs sum of the initial offset
    //Joint Scale And Offset Estimation

    //forward project
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
          //yEstimate->d[i_theta][i_r][i_t]=0;
          Real_t ttmp = sinogram->counts->getValue(i_theta, i_r, i_t) - errorSinogram->getValue(i_theta, i_r, i_t) - m_Mu->d[i_theta];
          yEstimate->setValue(ttmp, i_theta, i_r, i_t);
          yEstimate->divideByValue(m_I_0->d[i_theta], i_theta, i_r, i_t);
        }
      }
    }

    START_TIMER;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      Real_t a = 0;
      Real_t b = 0;
      Real_t c = 0;
      Real_t d = 0;
//  DATA_TYPE e = 0;
      Real_t numerator_sum = 0;
      Real_t denominator_sum = 0;
//  DATA_TYPE temp = 0.0;

      //compute the parameters of the quadratic for each view
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
          size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
          size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
          size_t yest_idx = yEstimate->calcIndex(i_theta, i_r, i_t);

          numerator_sum += (sinogram->counts->d[counts_idx] * m_Weight->d[weight_idx]);
          denominator_sum += (m_Weight->d[weight_idx]);

          a += (yEstimate->d[yest_idx] * m_Weight->d[weight_idx]);
          b += (yEstimate->d[yest_idx] * m_Weight->d[weight_idx] * sinogram->counts->d[counts_idx]);
          c += (sinogram->counts->d[counts_idx] * sinogram->counts->d[counts_idx] * m_Weight->d[weight_idx]);
          d += (yEstimate->d[yest_idx] * yEstimate->d[yest_idx] * m_Weight->d[weight_idx]);

        }
      }

      m_BkCost->setValue(numerator_sum, i_theta, 1); //yt*\lambda*1
      m_BkCost->setValue(b, i_theta, 0); //yt*\lambda*(Ax)
      m_CkCost->d[i_theta] = c; //yt*\lambda*y
      m_QkCost->setValue(denominator_sum, i_theta, 2);
      m_QkCost->setValue(a, i_theta, 1);
      m_QkCost->setValue(d, i_theta, 0);

      m_D1->d[i_theta] = numerator_sum / denominator_sum;
      m_D2->d[i_theta] = a / denominator_sum;

      a = 0;
      b = 0;
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
          size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
          size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
          size_t yest_idx = yEstimate->calcIndex(i_theta, i_r, i_t);

          a += ((yEstimate->d[yest_idx] - m_D2->d[i_theta]) * m_Weight->d[weight_idx] * yEstimate->d[yest_idx]);
          b -= ((sinogram->counts->d[counts_idx] - m_D1->d[i_theta]) * m_Weight->d[weight_idx] * yEstimate->d[yest_idx]);
        }
      }
      m_QuadraticParameters->setValue(a, i_theta, 0);
      m_QuadraticParameters->setValue(b, i_theta, 1);

#if 0
      temp = (m_QuadraticParameters->getValue(i_theta, 1) * m_QuadraticParameters->getValue(i_theta, 1)) / (4 * m_QuadraticParameters->getValue(i_theta, 0));

      if(temp > 0 && temp < high)
      {
        high = temp;
      } //high holds the maximum value for the rooting operation. beyond this value discriminants become negative. Basically high = min{b^2/4*a}
      else if(temp < 0 && temp > low)
      {
        low = temp;
      }
#endif
    }
    STOP_TIMER;
    PRINT_TIME("Joint Estimation Loops Time");

#ifdef COST_CALCULATE
    //compute cost
    /********************************************************************************************/
    Real_t sum = 0;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      sum += (m_QkCost->getValue(i_theta, 0) * m_I_0->d[i_theta] * m_I_0->d[i_theta]
          + 2 * m_QkCost->getValue(i_theta, 1) * m_I_0->d[i_theta] * m_Mu->d[i_theta]
          + m_Mu->d[i_theta] * m_Mu->d[i_theta] * m_QkCost->getValue(i_theta, 2)
          - 2 * (m_BkCost->getValue(i_theta, 0) * m_I_0->d[i_theta] + m_Mu->d[i_theta] * m_BkCost->getValue(i_theta, 1)) + m_CkCost->d[i_theta]); //evaluating the cost function
    }
    sum /= 2;
    printf("The value of the data match error prior to updating the I and mu =%lf\n", sum);

    /********************************************************************************************/

#endif //Cost calculate
    Real_t sum1 = 0;
    Real_t sum2 = 0;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      sum1 += (1.0 / (m_QkCost->getValue(i_theta, 0) - m_QkCost->getValue(i_theta, 1) * m_D2->d[i_theta]));
      sum2 += ((m_BkCost->getValue(i_theta, 0) - m_QkCost->getValue(i_theta, 1) * m_D1->d[i_theta])
          / (m_QkCost->getValue(i_theta, 0) - m_QkCost->getValue(i_theta, 1) * m_D2->d[i_theta]));
    }
    Real_t LagrangeMultiplier = (-sinogram->N_theta * m_TargetGain + sum2) / sum1;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {

      AverageMagI_k += fabs(m_I_0->d[i_theta]); //store the sum of the vector of gains

      Real_t NewI_k = (-1 * LagrangeMultiplier - m_QkCost->getValue(i_theta, 1) * m_D1->d[i_theta] + m_BkCost->getValue(i_theta, 0))
          / (m_QkCost->getValue(i_theta, 0) - m_QkCost->getValue(i_theta, 1) * m_D2->d[i_theta]);

      AverageI_kUpdate += fabs(NewI_k - m_I_0->d[i_theta]);

      m_I_0->d[i_theta] = NewI_k;
      //Postivity Constraint on the gains

      if(m_I_0->d[i_theta] < 0)
      {
        m_I_0->d[i_theta] *= 1;
      }
      AverageMagDelta_k += fabs(m_Mu->d[i_theta]);

      Real_t NewDelta_k = m_D1->d[i_theta] - m_D2->d[i_theta] * m_I_0->d[i_theta]; //some function of I_0[i_theta]
      AverageDelta_kUpdate += fabs(NewDelta_k - m_Mu->d[i_theta]);
      m_Mu->d[i_theta] = NewDelta_k;
      //Postivity Constraing on the offsets

      if(m_Mu->d[i_theta] < 0)
      {
        m_Mu->d[i_theta] *= 1;
      }
    }

#ifdef COST_CALCULATE
    /********************************************************************************************/
    //checking to see if the cost went down
    sum = 0;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      sum += (m_QkCost->getValue(i_theta, 0) * m_I_0->d[i_theta] * m_I_0->d[i_theta])
      + (2 * m_QkCost->getValue(i_theta, 1) * m_I_0->d[i_theta] * m_Mu->d[i_theta])
      + (m_Mu->d[i_theta] * m_Mu->d[i_theta] * m_QkCost->getValue(i_theta, 2))
      - (2 * (m_BkCost->getValue(i_theta, 0) * m_I_0->d[i_theta] + m_Mu->d[i_theta] * m_BkCost->getValue(i_theta, 1)) + m_CkCost->d[i_theta]); //evaluating the cost function
    }
    sum /= 2;

    printf("The value of the data match error after updating the I and mu =%lf\n", sum);
    /*****************************************************************************************************/
#endif //Cost calculate
    //Reproject to compute Error Sinogram for ICD
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
          size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
          size_t yest_idx = yEstimate->calcIndex(i_theta, i_r, i_t);
          size_t error_idx = errorSinogram->calcIndex(i_theta, i_r, i_t);

          errorSinogram->d[error_idx] = sinogram->counts->d[counts_idx] - m_Mu->d[i_theta] - (m_I_0->d[i_theta] * yEstimate->d[yest_idx]);
        }
      }
    }


    if(getVeryVerbose())
    {
      ss.str("");
      ss << "Lagrange Multiplier = " << LagrangeMultiplier;

      std::cout << "Tilt\tGains\tOffsets" << std::endl;
      for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
      {
        std::cout << i_theta << "\t" << m_I_0->d[i_theta] << "\t" << m_Mu->d[i_theta] << std::endl;
      }
      //std::cout << "Ratio of change in I_k " << AverageI_kUpdate / AverageMagI_k << std::endl;
      //std::cout << "Ratio of change in Delta_k " << AverageDelta_kUpdate / AverageMagDelta_k << std::endl;
    }
  }
#else 
  {
      Real_t num_sum = 0;
      Real_t den_sum = 0;
      Real_t alpha = 0;
	  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	  {
      for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
        {
		  if(m_Selector->getValue(i_theta, i_r, i_t) == 1)
		  {
          num_sum += (errorSinogram->getValue(i_theta, i_r, i_t) * m_Weight->getValue(i_theta, i_r, i_t));
          den_sum += m_Weight->getValue(i_theta, i_r, i_t);
		  }
        }
      }
	  }
      alpha = num_sum / den_sum;

	  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	  {	  
	     for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
         {
           for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
           {
          errorSinogram->deleteFromValue(alpha, i_theta, i_r, i_t);
           }
         }
      m_Mu->d[i_theta] += alpha;
	  if(getVeryVerbose())
	  {
	    std::cout << "Theta: " << i_theta << " Mu: " << m_Mu->d[i_theta] << std::endl;
	  }
	  }
     		  
  } 
#endif //BF_RECON
  //return 0;
}


// -----------------------------------------------------------------------------
// Updating the Weights for Noise Model
// -----------------------------------------------------------------------------
void HAADFForwardModel::updateWeights(SinogramPtr sinogram, RealVolumeType::Pointer ErrorSino)
{
  Real_t AverageVarUpdate = 0; //absolute sum of the gain updates
  Real_t AverageMagVar = 0; //absolute sum of the initial gains
  Real_t sum = 0;

  for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
  {
    sum = 0;
    //Factoring out the variance parameter from the Weight matrix
    for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
      {
        size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
        //     size_t yest_idx = yEstimate->calcIndex(i_theta, i_r, i_t);
        //    size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
#ifndef IDENTITY_NOISE_MODEL
        size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
  /*      if(sinogram->counts->d[counts_idx] != 0)
        {
          m_Weight->d[weight_idx] = 1.0 / sinogram->counts->d[counts_idx];
        }
        else
        {
          m_Weight->d[weight_idx] = 1.0;
        }*/
#ifdef BF_RECON //Override old weights
		  m_Weight->d[weight_idx] = BF_MAX/exp(sinogram->counts->d[counts_idx]);  
#endif //BF_RECON
		  
#else
        m_Weight->d[weight_idx] = 1.0;
#endif//Identity noise Model
      }

    }

    for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
      {
        size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
        size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
#ifdef BF_RECON
		if(m_Selector->d[weight_idx] == 1)  
#endif //BF_RECON
        sum += (ErrorSino->d[error_idx] * ErrorSino->d[error_idx] * m_Weight->d[weight_idx]); //Changed to only account for the counts
      }
    }
    sum /= (sinogram->N_r * sinogram->N_t);

    AverageMagVar += fabs(m_Alpha->d[i_theta]);
    AverageVarUpdate += fabs(sum - m_Alpha->d[i_theta]);
    m_Alpha->d[i_theta] = sum;
    //Update the weight for ICD updates
    for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
      {

        size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
#ifndef IDENTITY_NOISE_MODEL
        size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
  /*      if(m_Alpha->d[i_theta] != 0 && sinogram->counts->d[counts_idx] != 0)
        {
          m_Weight->d[weight_idx] = 1.0 / (sinogram->counts->d[counts_idx] * m_Alpha->d[i_theta]);
        }
        else
        {
          m_Weight->d[weight_idx] = 1.0;
        }*/
		
#ifdef BF_RECON
		  m_Weight->d[weight_idx] = (BF_MAX/exp(sinogram->counts->d[counts_idx]))/m_Alpha->d[i_theta];  
#endif //BF_RECON
		  
#else
        m_Weight->d[weight_idx] = 1.0 / m_Alpha->d[i_theta];
#endif //IDENTITY_NOISE_MODEL endif
	
      
	  }
    }

  }

	// This block averages the alphas so we have a single constant of proportionality for BF case 
	sum=0;		  
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	{
		sum+=m_Alpha->d[i_theta];
	}
	sum/=sinogram->N_theta;
	
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	{
		m_Alpha->d[i_theta]=sum;
	}
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	    for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
			for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
			{
				size_t weight_idx = m_Weight->calcIndex(i_theta, i_r, i_t);
#ifndef IDENTITY_NOISE_MODEL
				m_Weight->d[weight_idx] = (BF_MAX/exp(sinogram->counts->d[weight_idx]))/m_Alpha->d[i_theta];
#else
				m_Weight->d[weight_idx] = 1.0/m_Alpha->d[i_theta];
#endif //IDENTITY_NOISE_MODEL
			}
		
	
  if(getVeryVerbose())
  {
    std::cout << "Noise Model Weights:" << std::endl;
    std::cout << "Tilt\tWeight" << std::endl;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      std::cout << i_theta << "\t" << m_Alpha->d[i_theta] << std::endl;
    }
    std::cout << "Ratio of change in Variance " << AverageVarUpdate / AverageMagVar << std::endl;
  }

  notify("Update Weights Complete", 0, Observable::UpdateProgressMessage);
}

#ifdef BF_RECON

// -----------------------------------------------------------------------------
// Updating the boolean selector based on the error and weights
// -----------------------------------------------------------------------------
void HAADFForwardModel::updateSelector(SinogramPtr sinogram,
					RealVolumeType::Pointer ErrorSino)
{
for (uint16_t i_theta = 0; i_theta < sinogram->N_theta;i_theta++)
	for (uint16_t i_r = 0; i_r < sinogram->N_r; i_r++)
		for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
		{
			size_t idx = m_Weight->calcIndex(i_theta, i_r, i_t);
			if(ErrorSino->d[idx] * ErrorSino->d[idx] * m_Weight->d[idx] < m_BraggThreshold*m_BraggThreshold)
				m_Selector->d[idx] = 1;
			else 
				m_Selector->d[idx] = 0;
		}

}
#endif
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::writeNuisanceParameters(SinogramPtr sinogram)
{
  NuisanceParamWriter::Pointer nuisanceBinWriter = NuisanceParamWriter::New();
  nuisanceBinWriter->setNtheta(sinogram->N_theta);

  if(m_AdvParams->JOINT_ESTIMATION)
  {
    nuisanceBinWriter->setFileName(m_TomoInputs->gainsOutputFile);
    nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_I_O);
    nuisanceBinWriter->setData(m_I_0);
    nuisanceBinWriter->execute();
    if(nuisanceBinWriter->getErrorCondition() < 0)
    {
      setErrorCondition(-1);
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
    }

    nuisanceBinWriter->setFileName(m_TomoInputs->offsetsOutputFile);
    nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_mu);
    nuisanceBinWriter->setData(m_Mu);
    nuisanceBinWriter->execute();
    if(nuisanceBinWriter->getErrorCondition() < 0)
    {
      setErrorCondition(-1);
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
    }
  }

  if(m_AdvParams->NOISE_ESTIMATION)
  {
    nuisanceBinWriter->setFileName(m_TomoInputs->varianceOutputFile);
    nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_alpha);
    nuisanceBinWriter->setData(m_Alpha);
    nuisanceBinWriter->execute();
    if(nuisanceBinWriter->getErrorCondition() < 0)
    {
      setErrorCondition(-1);
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
    }
  } //Noise Model

  if(getVerbose())
  {
    std::cout << "Tilt\tFinal Gains\tFinal Offsets\tFinal Variances" << std::endl;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {

      if(m_AdvParams->NOISE_ESTIMATION)
      {
        std::cout << i_theta << "\t" << m_I_0->d[i_theta] << "\t" << m_Mu->d[i_theta] << "\t" << m_Alpha->d[i_theta] << std::endl;
      }
      else
      {
        std::cout << i_theta << "\t" << m_I_0->d[i_theta] << "\t" << m_Mu->d[i_theta] << std::endl;
      }
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::writeSinogramFile(SinogramPtr sinogram, RealVolumeType::Pointer finalSinogram)
{
  // Write the Sinogram out to a file
  SinogramBinWriter::Pointer sinogramWriter = SinogramBinWriter::New();
  sinogramWriter->setSinogram(sinogram);
  sinogramWriter->setTomoInputs(m_TomoInputs);
  sinogramWriter->setAdvParams(m_AdvParams);
  sinogramWriter->setObservers(getObservers());
  sinogramWriter->setI_0(m_I_0);
  sinogramWriter->setMu(m_Mu);
  sinogramWriter->setData(finalSinogram);
  sinogramWriter->execute();
  if(sinogramWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    notify(sinogramWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::createNuisanceParameters(SinogramPtr sinogram)
{
  int err = 0;
  err |= createInitialGainsData(sinogram);
  err |= createInitialOffsetsData(sinogram);
  err |= createInitialVariancesData(sinogram);
  return err;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::createInitialGainsData(SinogramPtr sinogram)
{
  std::stringstream ss;
  /* ********************* Initialize the Gains Array **************************/
  size_t gains_dims[1] =
  { sinogram->N_theta };
  m_InitialGain = RealArrayType::New(gains_dims, "sinogram->InitialGain");
  if(m_TomoInputs->gainsInputFile.empty() == false)
  {
    // Read the initial Gains from a File
    NuisanceParamReader::Pointer gainsInitializer = NuisanceParamReader::New();
    gainsInitializer->setFileName(m_TomoInputs->gainsInputFile);
    gainsInitializer->setData(m_InitialGain);
    gainsInitializer->setSinogram(sinogram);
    gainsInitializer->setAdvParams(m_AdvParams);
    gainsInitializer->setTomoInputs(m_TomoInputs);
    //   gainsInitializer->setGeometry(geometry);
    gainsInitializer->setObservers(getObservers());
    gainsInitializer->execute();
    if(gainsInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Gains from Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(gainsInitializer->getErrorCondition());
      return -1;
    }
  }
  else
  {
    // Set the values to the target gain value set by the user
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      m_InitialGain->d[i_theta] = m_TargetGain;
    }
  }
  /********************REMOVE************************/
  ss << "HARD WIRED TARGET GAIN" << std::endl;
  ss << "Target Gain: " << m_TargetGain << std::endl;
  /*************************************************/

  if(getVeryVerbose())
  {
    std::cout << ss.str() << std::endl;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::createInitialOffsetsData(SinogramPtr sinogram)
{
  std::stringstream ss;
  /* ********************* Initialize the Offsets Array **************************/
  size_t offsets_dims[1] =
  { sinogram->N_theta };
  m_InitialOffset = RealArrayType::New(offsets_dims, "sinogram->InitialOffset");
  if(m_TomoInputs->offsetsInputFile.empty() == false)
  {
    // Read the initial offsets from a File
    NuisanceParamReader::Pointer offsetsInitializer = NuisanceParamReader::New();
    offsetsInitializer->setFileName(m_TomoInputs->offsetsInputFile);
    offsetsInitializer->setData(m_InitialOffset);
    offsetsInitializer->setSinogram(sinogram);
    offsetsInitializer->setAdvParams(m_AdvParams);
    offsetsInitializer->setTomoInputs(m_TomoInputs);
//    offsetsInitializer->setGeometry(geometry);
    offsetsInitializer->setObservers(getObservers());
    offsetsInitializer->execute();
    if(offsetsInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Offsets from Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(offsetsInitializer->getErrorCondition());
      return -1;
    }
  }
  else if(m_UseDefaultOffset == true)
  {
    // Set the values to the default offset value set by the user
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      m_InitialOffset->d[i_theta] = m_DefaultOffset;
    }
  }
  else
  {
    // Compute the initial offset values from the data
    ComputeInitialOffsets::Pointer initializer = ComputeInitialOffsets::New();
    initializer->setSinogram(sinogram);
    initializer->setTomoInputs(m_TomoInputs);
    initializer->setAdvParams(m_AdvParams);
    initializer->setInitialOffset(m_InitialOffset);
    initializer->setObservers(getObservers());
    initializer->setVerbose(getVerbose());
    initializer->setVeryVerbose(getVeryVerbose());
    initializer->execute();
    if(initializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Offsets Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(initializer->getErrorCondition());
      return -1;
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADFForwardModel::createInitialVariancesData(SinogramPtr sinogram)
{

  /* ********************* Initialize the Variances Array **************************/
  size_t variance_dims[1] =
  { sinogram->N_theta };
  m_InitialVariance = RealArrayType::New(variance_dims, "sinogram->InitialVariance");
  if(m_TomoInputs->varianceInputFile.empty() == false)
  {
    // Read the initial variances from a File
    NuisanceParamReader::Pointer variancesInitializer = NuisanceParamReader::New();
    variancesInitializer->setFileName(m_TomoInputs->varianceInputFile);
    variancesInitializer->setData(m_InitialVariance);
    variancesInitializer->setSinogram(sinogram);
    variancesInitializer->setTomoInputs(m_TomoInputs);
    //   variancesInitializer->setGeometry(geometry);
    variancesInitializer->setAdvParams(m_AdvParams);
    variancesInitializer->setObservers(getObservers());
    variancesInitializer->setVerbose(getVerbose());
    variancesInitializer->setVeryVerbose(getVeryVerbose());
    variancesInitializer->execute();
    if(variancesInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Variances from Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(variancesInitializer->getErrorCondition());
      return -1;
    }
  }
  else
  {
    std::stringstream ss;
    for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
    {
      m_InitialVariance->d[i_theta] = m_DefaultVariance;
      ss << "Tilt: " << i_theta << "  Variance: " << m_InitialVariance->d[i_theta] << std::endl;
    }
    if(getVeryVerbose())
    {
      std::cout << ss.str() << std::endl;
    }
  }

  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Real_t HAADFForwardModel::forwardCost(SinogramPtr sinogram,RealVolumeType::Pointer ErrorSino)
{
	Real_t cost = 0,temp=0;
	Real_t errSinoValue = 0.0;
	
	
	//Data Mismatch Error
	
	for (int16_t i = 0; i < sinogram->N_theta; i++)
	{
		for (int16_t j = 0; j < sinogram->N_r; j++)
		{
			for (int16_t k = 0; k < sinogram->N_t; k++)
			{
				errSinoValue = ErrorSino->getValue(i, j, k);
				if(m_Selector->getValue(i,j,k) == 1)
				    cost += (errSinoValue * errSinoValue * m_Weight->getValue(i, j, k));
				else
					cost += m_BraggThreshold*m_BraggThreshold;	
			}
		}
	}
	
	cost /= 2;
	
	//Noise Error
	if(m_AdvParams->NOISE_ESTIMATION)
	{
		temp = 0;
		for (int16_t i = 0; i < sinogram->N_theta; i++)
		{
			for (int16_t j = 0; j < sinogram->N_r; j++)
			{
				for (int16_t k = 0; k < sinogram->N_t; k++)
				{
					if(m_Weight->getValue(i, j, k) != 0) temp += log(2 * M_PI * (1.0 / m_Weight->getValue(i, j, k)));
				}
			}
		}
		temp /= 2;
		cost += temp;
	} //NOISE_MODEL
	
	return cost;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::computeTheta(size_t Index,
					 std::vector<HAADFAMatrixCol::Pointer> &TempCol,
					 int32_t xzSliceIdx,
					 std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
					 RealVolumeType::Pointer ErrorSino,
					 SinogramPtr sinogram,
					 Int32ArrayType::Pointer Thetas)
{
	
	Thetas->d[0]=0;
	Thetas->d[1]=0;
	
	for (uint32_t q = 0; q < TempCol[Index]->count; q++)
	{
		uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (sinogram->N_r)));
		uint16_t i_r = (TempCol[Index]->index[q] % (sinogram->N_r));
		Real_t kConst0 = m_I_0->d[i_theta] * (TempCol[Index]->values[q]);
		uint16_t VoxelLineAccessCounter = 0;
		uint32_t vlrCount = VoxelLineResponse[xzSliceIdx]->index[0] + VoxelLineResponse[xzSliceIdx]->count;
		for (uint32_t i_t = VoxelLineResponse[xzSliceIdx]->index[0]; i_t < vlrCount; i_t++)
		{

			size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
			if(m_Selector->d[error_idx] == 1)
			{
			Real_t ProjectionEntry = kConst0 * VoxelLineResponse[xzSliceIdx]->values[VoxelLineAccessCounter];
			if(getBF_Flag() == false)
			{
				Thetas->d[1] += (ProjectionEntry * ProjectionEntry * m_Weight->d[error_idx]);
				Thetas->d[0] += (ErrorSino->d[error_idx] * ProjectionEntry * m_Weight->d[error_idx]);
			}
			else
			{
				ProjectionEntry *= getBFSinogram()->counts->d[error_idx];
				Thetas->d[1] += (ProjectionEntry * ProjectionEntry * m_Weight->d[error_idx]);
				Thetas->d[0] += (ErrorSino->d[error_idx] * ProjectionEntry * m_Weight->d[error_idx]);
			}
			VoxelLineAccessCounter++;
			}
		}
	}	
	Thetas->d[0]*=-1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::updateErrorSinogram(Real_t ChangeInVoxelValue,
						   size_t Index,
						   std::vector<HAADFAMatrixCol::Pointer> &TempCol,
						   int32_t xzSliceIdx,
						   std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
						   RealVolumeType::Pointer ErrorSino,
						   SinogramPtr sinogram)
{
	Real_t kConst2 = 0.0;
	//Update the ErrorSinogram
	for (uint32_t q = 0; q < TempCol[Index]->count; q++)
	{
		uint16_t i_theta = floor(static_cast<float>(TempCol[Index]->index[q] / (sinogram->N_r)));
		uint16_t i_r = (TempCol[Index]->index[q] % (sinogram->N_r));
		uint16_t VoxelLineAccessCounter = 0;
		for (uint32_t i_t = VoxelLineResponse[xzSliceIdx]->index[0]; i_t < VoxelLineResponse[xzSliceIdx]->index[0] + VoxelLineResponse[xzSliceIdx]->count; i_t++)
		{
			size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
			kConst2 = (m_I_0->d[i_theta]
                       * (TempCol[Index]->values[q] * VoxelLineResponse[xzSliceIdx]->values[VoxelLineAccessCounter] * (ChangeInVoxelValue)));
			if(getBF_Flag() == false)
			{
				ErrorSino->d[error_idx] -= kConst2;
			}
			else
			{
				
				ErrorSino->d[error_idx] -= (getBFSinogram()->counts->d[error_idx] * kConst2);
			}
			VoxelLineAccessCounter++;
#ifdef BF_RECON
        // if(ErrorSino->d[error_idx]*ErrorSino->d[error_idx]*m_Weight->d[error_idx] < m_BraggThreshold*m_BraggThreshold)
		   if(fabs(ErrorSino->d[error_idx]*sqrt(m_Weight->d[error_idx])) < m_BraggThreshold)
			 m_Selector->d[error_idx] = 1;
		   else {
			m_Selector->d[error_idx] = 0;	
			}
#endif
		}
	}
}


#ifdef BF_RECON
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
/*void HAADFForwardModel::setUpBraggThreshold(Real_t Threshold)
{
	m_BraggThreshold = Threshold;
}*/

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::processRawCounts(SinogramPtr sinogram)
{
    Real_t mean=0;
    for (int16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++) //slice index
    {
        for (int16_t i_r = 0; i_r < sinogram->N_r; i_r++)
        {
            for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
            {
                size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
                sinogram->counts->d[counts_idx] += BF_OFFSET;
                sinogram->counts->d[counts_idx] = -log(sinogram->counts->d[counts_idx]/BF_MAX);
				
               // if(sinogram->counts->d[counts_idx] < 0 ) //Clip the log data to be positive
               //     sinogram->counts->d[counts_idx] = 0;
				
                mean+=sinogram->counts->d[counts_idx];
            }
        }
    }
    mean/=(sinogram->N_theta*sinogram->N_r*sinogram->N_t);
    std::cout<<"Mean log value ="<<mean<<std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADFForwardModel::printRatioSelected(SinogramPtr sinogram)
{
	Real_t sum=0;
    for (int16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++) //slice index
    {
        for (int16_t i_r = 0; i_r < sinogram->N_r; i_r++)
        {
            for (uint16_t i_t = 0; i_t < sinogram->N_t; i_t++)
            {
			    size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
				sum+=m_Selector->d[counts_idx];
			}
		}
	}
	
	std::cout<<"Ratio of singoram entries used="<<sum/(sinogram->N_theta*sinogram->N_r*sinogram->N_t)<<std::endl;
}

void HAADFForwardModel::writeSelectorMrc(SinogramPtr sinogram,GeometryPtr geometry,RealVolumeType::Pointer ErrorSino)
{
	const std::string mrcFile="Selector.mrc";
	geometry->N_x = sinogram->N_r;
	geometry->N_y = sinogram->N_t;
	geometry->N_z = sinogram->N_theta;
	geometry->Object = sinogram->counts;
	for(uint32_t i_theta = 0; i_theta <  geometry->N_z; i_theta++)
	{
	   for(uint32_t i_r = 0; i_r < geometry->N_x; i_r++)
		  for(uint32_t i_t = 0; i_t < geometry->N_y; i_t++)
		  {
			size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
			Real_t value = exp(-sinogram->counts->d[counts_idx]+ErrorSino->d[counts_idx]);//m_Selector->d[counts_idx];
			value -= BF_OFFSET;
			counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);  
			geometry->Object->d[counts_idx] = value;			  
		  }
	}
	uint16_t cropStart=0;
	uint16_t cropEnd = geometry->N_x;
	/* Write the output to the MRC File */
	std::stringstream ss;
	ss.str("");
	ss << "Writing selector MRC file to '" << mrcFile << "'";
	notify(ss.str(), 0, Observable::UpdateProgressMessage);
	
	MRCWriter::Pointer mrcWriter = MRCWriter::New();
	mrcWriter->setOutputFile(mrcFile);
	mrcWriter->setGeometry(geometry);
	mrcWriter->setAdvParams(m_AdvParams);
	mrcWriter->setXDims(cropStart, cropEnd);
	mrcWriter->setYDims(0, geometry->N_y);
	mrcWriter->setZDims(0, geometry->N_z);
	mrcWriter->setObservers(getObservers());
	mrcWriter->execute();
	if(mrcWriter->getErrorCondition() < 0)
	{
		ss.str("");
		ss << "Error writing MRC file\n    '" << mrcFile << "'" << std::endl;
		setErrorCondition(mrcWriter->getErrorCondition());
		notify(ss.str(), 0, Observable::UpdateErrorMessage);
	}
}

Real_t HAADFForwardModel::estimateBraggThreshold(SinogramPtr sinogram, RealVolumeType::Pointer ErrorSino,Real_t percentage)
{
	//m_Alpha->d[i_theta]
	//m_Weight->d[error_idx]
	Real_t EstBraggThresh;
	uint32_t NumElts = sinogram->N_theta*sinogram->N_r*sinogram->N_t;
	uint32_t NumEltsReject = percentage*NumElts;//Find the Error*Weight 
	//corresponding to this order statistic
	RealArrayType::Pointer Ratio;
	size_t dims[1]={NumElts};
	Ratio = RealArrayType::New(dims, "Ratio of error to noise variance");
	uint32_t counts=0;
	for(uint32_t i_theta = 0; i_theta <  sinogram->N_theta; i_theta++)
		for(uint32_t i_r = 0; i_r < sinogram->N_r; i_r++)
			for(uint32_t i_t = 0; i_t < sinogram->N_t; i_t++)
			{
				size_t counts_idx = sinogram->counts->calcIndex(i_theta, i_r, i_t);
				Ratio->d[counts]=ErrorSino->d[counts_idx]*m_Weight->d[counts_idx];
				Ratio->d[counts]*=ErrorSino->d[counts_idx];
				counts++;
			}
	std::cout<<"Num Elts"<<counts<<std::endl;
	std::cout<<"Num Elts to reject ="<<NumEltsReject<<std::endl;
	
	EstBraggThresh =sqrt(RandomizedSelect(Ratio,0, counts-1, NumElts - NumEltsReject));
	std::cout<<"Bragg Thresh estimated using Randomized select"<<EstBraggThresh<<std::endl;	
	
	/*uint32_t max_index=0;
	for(uint32_t j =0; j < NumEltsReject;j++)
	{
	Real_t max=-INFINITY;	
	for(uint32_t i =j ;i < NumElts; i++)
	{
		if (Ratio->d[i] > max) {
			max =Ratio->d[i];
			max_index = i;
		}
	}
	Real_t temp = Ratio->d[j];
	Ratio->d[j] = max;
	Ratio->d[max_index]=temp;	
	}
	EstBraggThresh = sqrt(Ratio->d[NumEltsReject-1]);
	std::cout<<"Bragg Thresh estimated using Insertion sort = "<<EstBraggThresh<<std::endl;	
	*/
	return EstBraggThresh;
}

//Radomized select based on code from : http://stackoverflow.com/questions/5847273/order-statistic-implmentation


Real_t HAADFForwardModel::RandomizedSelect(RealArrayType::Pointer A,uint32_t p, uint32_t r,uint32_t i)
{
	if (p == r)
    {
		return A->d[p];
    }
	uint32_t q = RandomizedPartition(A, p, r);
	uint32_t k = q - p + 1;
	if (i == k)
    {
		return A->d[q];
    }
	else if  (i < k)
    {
		return RandomizedSelect(A, p, q-1, i) ;
    }
	else return RandomizedSelect(A, q+1, r, i - k);
}
uint32_t HAADFForwardModel::Partition(RealArrayType::Pointer A,uint32_t p,uint32_t r)
{
	Real_t x=A->d[r],temp;
	uint32_t i=p-1,j;
	for(j=p;j<r;j++)
    {
		if(A->d[j]<=x)
        {
			i++;
			temp=A->d[i];
			A->d[i]=A->d[j];
			A->d[j]=temp;
        }
    }
	temp=A->d[i+1];
	A->d[i+1]=A->d[r];
	A->d[r]=temp;
	return i+1;
	
}
uint32_t HAADFForwardModel::RandomizedPartition(RealArrayType::Pointer A,uint32_t p,uint32_t r)
{
	Real_t temp;
	uint32_t j = p + rand()%(r-p+1);
	temp = A->d[r];
	A->d[r] = A->d[j];
	A->d[j] = temp;	
	return Partition(A, p, r);
}

#endif //BF Recon