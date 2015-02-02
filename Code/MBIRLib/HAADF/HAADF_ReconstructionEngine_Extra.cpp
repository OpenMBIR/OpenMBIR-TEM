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

#include "HAADF_ReconstructionEngine.h"
#include "MBIRLib/HAADF/HAADFConstants.h"

// Read the Input data from the supplied data file
// We are scoping here so the various readers are automatically cleaned up before
// the code goes any farther
int HAADF_ReconstructionEngine::readInputData()
{
  TomoFilter::Pointer dataReader = TomoFilter::NullPointer();
  std::string extension = MXAFileInfo::extension(m_TomoInputs->sinoFile);
  if(extension.compare("bin") == 0)
  {
    dataReader = RawSinogramInitializer::NewTomoFilter();
  }
  else
    // if(extension.compare("mrc") == 0 || extension.compare("ali") == 0)
  {
    // We are going to assume that the user has selected a valid MRC file to get this far so we are just going to try
    // to read the file as an MRC file which may really cause issues but there does not seem to be any standard file
    // extensions for MRC files.
    dataReader = MRCSinogramInitializer::NewTomoFilter();
  }

  //  {
  //    setErrorCondition(-1);
  //    notify("A supported file reader for the input file was not found.", 100, Observable::UpdateErrorMessage);
  //    return -1;
  //  }
  dataReader->setTomoInputs(m_TomoInputs);
  dataReader->setSinogram(m_Sinogram);
  dataReader->setAdvParams(m_AdvParams);
  dataReader->setObservers(getObservers());
  dataReader->setVerbose(getVerbose());
  dataReader->setVeryVerbose(getVeryVerbose());
  dataReader->execute();
  if(dataReader->getErrorCondition() < 0)
  {
    notify("Error reading Input Sinogram Data file", 100, Observable::UpdateErrorMessage);
    setErrorCondition(dataReader->getErrorCondition());
    return -1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ReconstructionEngine::initializeBrightFieldData()
{
  std::stringstream ss;
  uint8_t Flag1 = 1;
  uint8_t Flag2 = 1;
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
      notify("A supported file reader for the Bright Field file was not found.", 100, Observable::UpdateErrorMessage);
      return -1;
    }
    dataReader->setTomoInputs(m_BFTomoInputs);
    dataReader->setSinogram(m_BFSinogram);
    dataReader->setAdvParams(m_AdvParams);
    dataReader->setObservers(getObservers());
    dataReader->execute();
    if(dataReader->getErrorCondition() < 0)
    {
      Flag1 = 0;
      //Try one more time - if the sizes are not matches just read in the
      //full BF sinogram
      m_BFTomoInputs->xStart = 0;
      m_BFTomoInputs->xEnd = m_Sinogram->N_r - 1;
      m_BFTomoInputs->yStart = 0;
      m_BFTomoInputs->yEnd = m_Sinogram->N_t - 1;
      m_BFTomoInputs->zStart = 0;
      m_BFTomoInputs->zEnd = m_Sinogram->N_theta - 1;

      dataReader->setTomoInputs(m_BFTomoInputs);
      dataReader->setSinogram(m_BFSinogram);
      dataReader->setAdvParams(m_AdvParams);
      dataReader->setObservers(getObservers());
      dataReader->execute();
      if(dataReader->getErrorCondition() < 0)
      {
        Flag2 = 0;
        //  notify("Error reading Input Sinogram Data file", 100, Observable::UpdateProgressValueAndMessage);
        //  setErrorCondition(dataReader->getErrorCondition());
        //  return -1;
      }
      if(Flag1 == 0 && Flag2 == 0)
      {
        notify("Error reading Input Sinogram Data file", 100, Observable::UpdateErrorMessage);
        setErrorCondition(dataReader->getErrorCondition());
        return -1;
      }
    }

    if(m_BFSinogram->N_r != m_Sinogram->N_r || m_BFSinogram->N_t != m_Sinogram->N_t || m_BFSinogram->N_theta != m_Sinogram->N_theta)
    {
      notify("The two file sizes are not matched", 100, Observable::UpdateErrorMessage);
      setErrorCondition(dataReader->getErrorCondition());
      return -1;
    }

    //Normalize the HAADF image
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      //slice index
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          //1000 is for Marc De Graef data which needed to multiplied
          //  Real_t ttmp = (m_BFSinogram->counts->getValue(i_theta, i_r, i_t) * 1000);
          //  m_Sinogram->counts->divideByValue(ttmp, i_theta, i_r, i_t);
          //100 is for Marc De Graef data which needed to multiplied
          //  m_BFSinogram->counts->multiplyByValue(100, i_theta, i_r, i_t);
          m_BFSinogram->counts->setValue(m_BFSinogram->counts->getValue(i_theta, i_r, i_t) + BF_OFFSET, i_theta, i_r, i_t);
        }
      }
    }
    m_ForwardModel->setBF_Flag(true);
    notify("BF initialization complete", 0, Observable::UpdateProgressMessage);
  }
  else
  {
    m_ForwardModel->setBF_Flag(false);
  }
  return 0;
}


// -----------------------------------------------------------------------------
// Initialize the Geometry data from a rough reconstruction
// -----------------------------------------------------------------------------
int HAADF_ReconstructionEngine::initializeRoughReconstructionData()
{
  InitialReconstructionInitializer::Pointer geomInitializer = InitialReconstructionInitializer::NullPointer();
  std::string extension = MXAFileInfo::extension(m_TomoInputs->initialReconFile);
  if (m_TomoInputs->initialReconFile.empty() == true)
  {
    // This will just initialize all the values to Zero (0) or a DefaultValue Set by user
    geomInitializer = InitialReconstructionInitializer::New();
  }
  else if (extension.compare("bin") == 0 )
  {
    // This will read the values from a binary file
    geomInitializer = InitialReconstructionBinReader::NewInitialReconstructionInitializer();
  }
  else if (extension.compare("mrc") == 0)
  {
    notify("We are not dealing with mrc volume files.", 0, Observable::UpdateErrorMessage);
    return -1;
  }
  else
  {
    notify("Could not find a compatible reader for the initial reconstruction data file. The program will now end.", 0, Observable::UpdateErrorMessage);
    return -1;
  }
  geomInitializer->setSinogram(m_Sinogram);
  geomInitializer->setTomoInputs(m_TomoInputs);
  geomInitializer->setGeometry(m_Geometry);
  geomInitializer->setAdvParams(m_AdvParams);
  geomInitializer->setObservers(getObservers());
  geomInitializer->setVerbose(getVerbose());
  geomInitializer->setVeryVerbose(getVeryVerbose());
  geomInitializer->execute();

  if(geomInitializer->getErrorCondition() < 0)
  {
    notify("Error reading Initial Reconstruction Data from File", 100, Observable::UpdateErrorMessage);
    setErrorCondition(geomInitializer->getErrorCondition());
    return -1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::initializeROIMask(UInt8Image_t::Pointer Mask)
{
  Real_t x = 0.0;
  Real_t z = 0.0;
  for (uint16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (uint16_t j = 0; j < m_Geometry->N_x; j++)
    {
      x = m_Geometry->x0 + ((Real_t)j + 0.5) * m_TomoInputs->delta_xz;
      z = m_Geometry->z0 + ((Real_t)i + 0.5) * m_TomoInputs->delta_xz;
      if(x >= -(m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && x <= (m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && z >= -m_TomoInputs->LengthZ / 2
          && z <= m_TomoInputs->LengthZ / 2)
      {
        Mask->setValue(1, i, j);
      }
      else
      {
        Mask->setValue(0, i, j);
      }
    }
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::computeOriginalXDims(uint16_t& cropStart, uint16_t& cropEnd)
{

  Real_t x;
  for (uint16_t j = 0; j < m_Geometry->N_x; j++)
  {
    x = m_Geometry->x0 + ((Real_t)j + 0.5) * m_TomoInputs->delta_xz;
    if(x >= -(m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && x <= (m_Sinogram->N_r * m_Sinogram->delta_r) / 2)
    {
      cropStart = j;
      break;
    }

  }

  for (int16_t j = m_Geometry->N_x - 1; j >= 0; j--)
  {
    x = m_Geometry->x0 + ((Real_t)j + 0.5) * m_TomoInputs->delta_xz;
    if(x >= -(m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && x <= (m_Sinogram->N_r * m_Sinogram->delta_r) / 2)
    {
      cropEnd = (uint16_t)j;
      break;
    }
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::gainAndOffsetInitialization()
{
  Real_t sum = 0;
  Real_t temp = 0;
  for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
  {
    // Gains
    m_ForwardModel->getI_0()->d[k] = m_ForwardModel->getInitialGain()->d[k];
    // Offsets
    m_ForwardModel->getMu()->d[k] = m_ForwardModel->getInitialOffset()->d[k];

    sum += m_ForwardModel->getI_0()->d[k];

  }
  sum /= m_Sinogram->N_theta;

  if (getVerbose()) { printf("The Arithmetic mean of the constraint is %lf\n", sum); }
  if(sum - m_ForwardModel->getTargetGain() > 1e-5)
  {
    if (getVerbose()) { printf("Arithmetic Mean Constraint not met..renormalizing\n");}
    temp = m_ForwardModel->getTargetGain() / sum;
    for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
    {
      m_ForwardModel->getI_0()->d[k] = m_ForwardModel->getInitialGain()->d[k] * temp;
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::initializeVolume(RealVolumeType::Pointer Y_Est, double value)
{
  for (uint16_t i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (uint16_t j = 0; j < m_Sinogram->N_r; j++)
    {
      for (uint16_t k = 0; k < m_Sinogram->N_t; k++)
      {
        Y_Est->setValue(value, i, j, k);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::storeVoxelResponse(RealVolumeType::Pointer H_t,
                                                    std::vector<AMatrixCol::Pointer>& VoxelLineResponse)
{
  Real_t ProfileThickness = 0.0;
  Real_t y = 0.0;
  Real_t t = 0.0;
  Real_t tmin;
  Real_t tmax;
  int16_t slice_index_min, slice_index_max;
  Real_t center_t, delta_t;
  int16_t index_delta_t;
  Real_t w3, w4;

  //Storing the response along t-direction for each voxel line
  notify("Storing the response along Y-direction for each voxel line", 0, Observable::UpdateProgressMessage);
  if(getVeryVerbose())
  {
    std::cout << "Voxel Response by Y Slice" << std::endl;
    std::cout << "Y\tProfile Thickness" << std::endl;
  }
  for (uint16_t i = 0; i < m_Geometry->N_y; i++)
  {
    y = ((Real_t)i + 0.5) * m_TomoInputs->delta_xy + m_Geometry->y0;
    t = y;
    tmin = (t - m_TomoInputs->delta_xy / 2) > m_Sinogram->T0 ? t - m_TomoInputs->delta_xy / 2 : m_Sinogram->T0;
    tmax = (t + m_TomoInputs->delta_xy / 2) <= m_Sinogram->TMax ? t + m_TomoInputs->delta_xy / 2 : m_Sinogram->TMax;

    slice_index_min = static_cast<uint16_t>(floor((tmin - m_Sinogram->T0) / m_Sinogram->delta_t));
    slice_index_max = static_cast<uint16_t>(floor((tmax - m_Sinogram->T0) / m_Sinogram->delta_t));

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
      center_t = ((Real_t)i_t + 0.5) * m_Sinogram->delta_t + m_Sinogram->T0;
      delta_t = fabs(center_t - t);
      index_delta_t = static_cast<uint16_t>(floor(delta_t / m_DetectorParameters->getOffsetT()));
      if(index_delta_t < m_AdvParams->DETECTOR_RESPONSE_BINS)
      {
        w3 = delta_t - (Real_t)(index_delta_t) *  m_DetectorParameters->getOffsetT();
        w4 = ((Real_t)index_delta_t + 1) *  m_DetectorParameters->getOffsetT() - delta_t;
        uint16_t ttmp = index_delta_t + 1 < m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : m_AdvParams->DETECTOR_RESPONSE_BINS - 1;
        ProfileThickness = (w4 / m_DetectorParameters->getOffsetT()) * H_t->getValue(0, 0, index_delta_t) + (w3 /  m_DetectorParameters->getOffsetT()) * H_t->getValue(0, 0, ttmp);
        //  ProfileThickness = (w4 / OffsetT) * detectorResponse->d[0][uint16_t(floor(m_Sinogram->N_theta/2))][index_delta_t]
        //  + (w3 / OffsetT) * detectorResponse->d[0][uint16_t(floor(m_Sinogram->N_theta/2))][index_delta_t + 1 <m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_t + 1 :m_AdvParams->DETECTOR_RESPONSE_BINS - 1];
      }
      else
      {
        ProfileThickness = 0;
      }

      if(ProfileThickness != 0) //Store the response of this slice
      {
        if(getVeryVerbose())
        {
          std::cout << i_t << "\t" << ProfileThickness << std::endl;
        }
        AMatrixCol::Pointer vlr = VoxelLineResponse[i];
        int32_t count = vlr->count;
        vlr->values[count] = ProfileThickness;
        vlr->index[count] = i_t;
        //        size_t dim0 = vlr->valuesPtr->getDims()[0];
        vlr->setCount(count + 1);
      }
    }
  }
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::calculateArithmeticMean()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ReconstructionEngine::jointEstimation(RealVolumeType::Pointer Weight,
                                                RealVolumeType::Pointer ErrorSino,
                                                RealVolumeType::Pointer Y_Est,
                                                CostData::Pointer cost)
{
  std::stringstream ss;
  std::string indent("  ");
  RealArrayType::Pointer I_0 = m_ForwardModel->getI_0();
  RealArrayType::Pointer mu = m_ForwardModel->getMu();
  RealArrayType::Pointer sigma = m_ForwardModel->getAlpha();
  //Estimate only offsets

    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      Real_t num_sum = 0;
      Real_t den_sum = 0;
      Real_t alpha = 0;
      Real_t temp =0;
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
	  size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
	  if(m_ForwardModel->getBraggSelector()->d[error_idx])
	    {
          num_sum += (ErrorSino->getValue(i_theta, i_r, i_t) * Weight->getValue(i_theta, i_r, i_t));
          den_sum += Weight->getValue(i_theta, i_r, i_t);
	    }
	  else
	    {
	      temp = m_ForwardModel->getBraggDelta()*m_ForwardModel->getBraggThreshold();
              temp /= abs(ErrorSino->getValue(i_theta,i_r,i_t));
              temp *= sqrt(Weight->getValue(i_theta,i_r,i_t)*sigma->d[i_theta]);
	      num_sum += temp*ErrorSino->getValue(i_theta,i_r,i_t);
	      den_sum += temp;
	    }
        }
      }
      alpha = num_sum / den_sum;

      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          ErrorSino->deleteFromValue(alpha, i_theta, i_r, i_t);
        }
      }
      mu->d[i_theta] += alpha;
      if(getVeryVerbose())
      {
        std::cout << "Theta: " << i_theta << " Mu: " << mu->d[i_theta] << std::endl;
      }
    }
    m_ForwardModel->computeBraggSelector(ErrorSino,Weight);
#ifdef COST_CALCULATE
    /*********************Cost Calculation*************************************/
    Real_t cost_value = computeCost(ErrorSino, Weight);
    std::cout << cost_value << std::endl;
    int increase = cost->addCostValue(cost_value);
    if (increase == 1)
    {
      std::cout << "Cost just increased after offset update!" << std::endl;
      //break;
      return -1;
    }
    cost->writeCostValue(cost_value);
    /**************************************************************************/
#endif

  return 0;
}
// -----------------------------------------------------------------------------
// Calculate Error Sinogram
// Also compute weights of the diagonal covariance matrix
// -----------------------------------------------------------------------------

void HAADF_ReconstructionEngine::calculateMeasurementWeight(RealVolumeType::Pointer Weight,
                                                            RealVolumeType::Pointer ErrorSino,
                                                            RealVolumeType::Pointer Y_Est)
{
  Real_t checksum = 0;
  START_TIMER;
  RealArrayType::Pointer alpha = m_ForwardModel->getAlpha();
  RealArrayType::Pointer mu = m_ForwardModel->getMu();

  for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++) //slice index
  {
    if (m_AdvParams->NOISE_MODEL)
    {
      alpha->d[i_theta] = m_ForwardModel->getInitialVariance()->d[i_theta]; //Initialize the refinement parameters from any previous run
    }//Noise model

    checksum = 0;
    for (int16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        size_t counts_idx = m_Sinogram->counts->calcIndex(i_theta, i_r, i_t);
        size_t weight_idx = Weight->calcIndex(i_theta, i_r, i_t);
        size_t yest_idx = Y_Est->calcIndex(i_theta, i_r, i_t);
        size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
        if(m_ForwardModel->getBF_Flag() == false)
        {
          ErrorSino->d[error_idx] = m_Sinogram->counts->d[counts_idx] - Y_Est->d[yest_idx] - mu->d[i_theta];
        }
        else
        {
          size_t bfcounts_idx = m_BFSinogram->counts->calcIndex(i_theta, i_r, i_t);

          ErrorSino->d[error_idx] = m_Sinogram->counts->d[counts_idx] - m_BFSinogram->counts->d[bfcounts_idx] * Y_Est->d[yest_idx]
                                    - mu->d[i_theta];
        }

#ifndef IDENTITY_NOISE_MODEL
        if(m_Sinogram->counts->d[counts_idx] != 0)
        {
          Weight->d[weight_idx] = 1.0 / m_Sinogram->counts->d[counts_idx];
        }
        else
        {
          Weight->d[weight_idx] = 1e-10; //Set the weight to some small number
          //TODO: Make this something resonable
        }
#else
        Weight->d[weight_idx] = 1.0;
#endif //IDENTITY_NOISE_MODEL endif
#ifdef FORWARD_PROJECT_MODE
        temp = Y_Est->d[i_theta][i_r][i_t] / I_0->d[i_theta];
        fwrite(&temp, sizeof(Real_t), 1, Fp6);
#endif
#ifdef DEBUG
        if(Weight->d[weight_idx] < 0)
        {
          //  std::cout << m_Sinogram->counts->d[counts_idx] << "    " << alpha->d[i_theta] << std::endl;
        }
#endif//Debug

        if (m_AdvParams->NOISE_MODEL)
        {
          Weight->d[weight_idx] /= alpha->d[i_theta];
        }// NOISE_MODEL


        checksum += Weight->d[weight_idx];
      }
    }
    if(getVerbose())
    {
      printf("Check sum of Diagonal Covariance Matrix= %lf\n", checksum);
    }
  }
  STOP_TIMER;
  std::string indent(" ");
  PRINT_TIME("Computing Weights");
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ReconstructionEngine::calculateCost(CostData::Pointer cost,
                                              RealVolumeType::Pointer Weight,
                                              RealVolumeType::Pointer ErrorSino)
{
  Real_t cost_value = computeCost(ErrorSino, Weight);
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
// Updating the Weights for Noise Model
// -----------------------------------------------------------------------------
void HAADF_ReconstructionEngine::updateWeights(RealVolumeType::Pointer Weight,
                                               RealVolumeType::Pointer ErrorSino)
{
  Real_t AverageVarUpdate = 0; //absolute sum of the gain updates
  Real_t AverageMagVar = 0; //absolute sum of the initial gains
  Real_t sum = 0;
  RealArrayType::Pointer alpha = m_ForwardModel->getAlpha();

  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
  {
    sum = 0;
    //Factoring out the variance parameter from the Weight matrix
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        size_t weight_idx = Weight->calcIndex(i_theta, i_r, i_t);
#ifndef IDENTITY_NOISE_MODEL
        size_t counts_idx = m_Sinogram->counts->calcIndex(i_theta, i_r, i_t);
        if(m_Sinogram->counts->d[counts_idx] != 0)
        {
          Weight->d[weight_idx] = 1.0 / m_Sinogram->counts->d[counts_idx];
        }
        else
        {
          Weight->d[weight_idx] = 1.0;
        }
#else
        Weight->d[weight_idx] = 1.0;
#endif//Identity noise Model
      }
    }

    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        size_t error_idx = ErrorSino->calcIndex(i_theta, i_r, i_t);
        if(m_ForwardModel->getBraggSelector()->d[error_idx])
	  sum += (ErrorSino->d[error_idx] * ErrorSino->d[error_idx] * Weight->d[error_idx]); //Changed to only account for the counts
	else
	  {
	    sum += m_ForwardModel->getBraggThreshold()*m_ForwardModel->getBraggDelta()*sqrt(Weight->d[error_idx])*alpha->d[i_theta];
	  }
      }
    }
    sum /= (m_Sinogram->N_r * m_Sinogram->N_t);

    AverageMagVar += fabs(alpha->d[i_theta]);
    AverageVarUpdate += fabs(sum - alpha->d[i_theta]);
    alpha->d[i_theta] = sum;
    //Update the weight for ICD updates
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {

        size_t weight_idx = Weight->calcIndex(i_theta, i_r, i_t);
#ifndef IDENTITY_NOISE_MODEL
        size_t counts_idx = m_Sinogram->counts->calcIndex(i_theta, i_r, i_t);
        if(alpha->d[i_theta] != 0 && m_Sinogram->counts->d[counts_idx] != 0)
        {
          Weight->d[weight_idx] = 1.0 / (m_Sinogram->counts->d[counts_idx] * alpha->d[i_theta]);
        }
        else
        {
          Weight->d[weight_idx] = 1.0;
        }
#else
        Weight->d[weight_idx] = 1.0 / alpha->d[i_theta];
#endif //IDENTITY_NOISE_MODEL endif
      }
    }

  }
  
  m_ForwardModel->computeBraggSelector(ErrorSino,Weight);

  if(getVeryVerbose())
  {
    std::cout << "Noise Model Weights:" << std::endl;
    std::cout << "Tilt\tWeight" << std::endl;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      std::cout << i_theta << "\t" << alpha->d[i_theta] << std::endl;
    }
    std::cout << "Ratio of change in Variance " << AverageVarUpdate / AverageMagVar << std::endl;
  }

  notify("Update Weights Complete", 0, Observable::UpdateProgressMessage);
}

#ifdef BF_RECON
void HAADF_ReconstructionEngine::processRawCounts()
{
  Real_t mean = 0;
  for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++) //slice index
  {
    for (int16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        size_t counts_idx = m_Sinogram->counts->calcIndex(i_theta, i_r, i_t);
        m_Sinogram->counts->d[counts_idx] += BF_OFFSET;
        m_Sinogram->counts->d[counts_idx] = -log(m_Sinogram->counts->d[counts_idx] / BF_MAX);

        if(m_Sinogram->counts->d[counts_idx] < 0 ) //Clip the log data to be positive
        { m_Sinogram->counts->d[counts_idx] = 0; }

        mean += m_Sinogram->counts->d[counts_idx];
      }
    }
  }
  mean /= (m_Sinogram->N_theta * m_Sinogram->N_r * m_Sinogram->N_t);
  std::cout << "Mean log value =" << mean << std::endl;
}

#endif //BF Recon




