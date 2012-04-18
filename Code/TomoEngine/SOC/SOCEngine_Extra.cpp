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

#include "SOCEngine.h"

// Read the Input data from the supplied data file
// We are scoping here so the various readers are automatically cleaned up before
// the code goes any farther
int SOCEngine::readInputData()
{
  TomoFilter::Pointer dataReader = TomoFilter::NullPointer();
  std::string extension = MXAFileInfo::extension(m_TomoInputs->sinoFile);
  if(extension.compare("mrc") == 0 || extension.compare("ali") == 0)
  {
    dataReader = MRCSinogramInitializer::NewTomoFilter();
  }
  else if(extension.compare("bin") == 0)
  {
    dataReader = RawSinogramInitializer::NewTomoFilter();
  }
  else
  {
    setErrorCondition(-1);
    notify("A supported file reader for the input file was not found.", 100, Observable::UpdateProgressValueAndMessage);
    return -1;
  }
  dataReader->setTomoInputs(m_TomoInputs);
  dataReader->setSinogram(m_Sinogram);
  dataReader->setObservers(getObservers());
  dataReader->execute();
  if(dataReader->getErrorCondition() < 0)
  {
    notify("Error reading Input Sinogram Data file", 100, Observable::UpdateProgressValueAndMessage);
    setErrorCondition(dataReader->getErrorCondition());
    return -1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int SOCEngine::initializeBrightFieldData()
{

  if (m_BFTomoInputs.get() != NULL && m_BFSinogram.get() != NULL && m_BFTomoInputs->sinoFile.empty()== false)
  {
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
    dataReader->setObservers(getObservers());
    dataReader->execute();
    if(dataReader->getErrorCondition() < 0)
    {
      notify("Error reading Input Sinogram Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(dataReader->getErrorCondition());
      return -1;
    }

    //Normalize the HAADF image
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    { //slice index
      for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      {
        for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
        {
          //1000 is for Marc De Graef data which needed to multiplied
          m_Sinogram->counts->d[i_theta][i_r][i_t] /= (m_BFSinogram->counts->d[i_theta][i_r][i_t] * 1000);
        }
      }
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int SOCEngine::createInitialGainsData()
{
  /* ********************* Initialize the Gains Array **************************/
  size_t gains_dims[1] =
  { m_Sinogram->N_theta };
  m_Sinogram->InitialGain = RealArrayType::New(gains_dims, "sinogram->InitialGain");
  if(m_TomoInputs->gainsInputFile.empty() == false)
  {
    // Read the initial Gains from a File
    NuisanceParamReader::Pointer gainsInitializer = NuisanceParamReader::New();
    gainsInitializer->setFileName(m_TomoInputs->gainsInputFile);
    gainsInitializer->setData(m_Sinogram->InitialGain);
    gainsInitializer->setSinogram(m_Sinogram);
    gainsInitializer->setTomoInputs(m_TomoInputs);
    gainsInitializer->setGeometry(m_Geometry);
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
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_Sinogram->InitialGain->d[i_theta] = m_TomoInputs->targetGain;
    }
  }
  /********************REMOVE************************/
  std::cout << "HARD WIRED TARGET GAIN" << std::endl;
  m_Sinogram->targetGain = m_TomoInputs->targetGain; //TARGET_GAIN;
  std::cout << "Target Gain: " << m_Sinogram->targetGain << std::endl;
  /*************************************************/

  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int SOCEngine::createInitialOffsetsData()
{
  /* ********************* Initialize the Offsets Array **************************/
  size_t offsets_dims[1] =
  { m_Sinogram->N_theta };
  m_Sinogram->InitialOffset = RealArrayType::New(offsets_dims, "sinogram->InitialOffset");
  if(m_TomoInputs->offsetsInputFile.empty() == false)
  {
    // Read the initial offsets from a File
    NuisanceParamReader::Pointer offsetsInitializer = NuisanceParamReader::New();
    offsetsInitializer->setFileName(m_TomoInputs->offsetsInputFile);
    offsetsInitializer->setData(m_Sinogram->InitialOffset);
    offsetsInitializer->setSinogram(m_Sinogram);
    offsetsInitializer->setTomoInputs(m_TomoInputs);
    offsetsInitializer->setGeometry(m_Geometry);
    offsetsInitializer->setObservers(getObservers());
    offsetsInitializer->execute();
    if(offsetsInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Offsets from Data file", 100, Observable::UpdateProgressValueAndMessage);
      setErrorCondition(offsetsInitializer->getErrorCondition());
      return -1;
    }
  }
  else if(m_TomoInputs->useDefaultOffset == true)
  {
    // Set the values to the default offset value set by the user
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_Sinogram->InitialOffset->d[i_theta] = m_TomoInputs->defaultOffset;
    }
  }
  else
  {
    // Compute the initial offset values from the data
    ComputeInitialOffsets::Pointer initializer = ComputeInitialOffsets::New();
    initializer->setSinogram(m_Sinogram);
    initializer->setTomoInputs(m_TomoInputs);
    initializer->setObservers(getObservers());
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
int SOCEngine::createInitialVariancesData()
{

  /* ********************* Initialize the Variances Array **************************/
  size_t variance_dims[1] =
  { m_Sinogram->N_theta };
  m_Sinogram->InitialVariance = RealArrayType::New(variance_dims, "sinogram->InitialVariance");
  if(m_TomoInputs->varianceInputFile.empty() == false)
  {
    // Read the initial variances from a File
    NuisanceParamReader::Pointer variancesInitializer = NuisanceParamReader::New();
    variancesInitializer->setFileName(m_TomoInputs->varianceInputFile);
    variancesInitializer->setData(m_Sinogram->InitialVariance);
    variancesInitializer->setSinogram(m_Sinogram);
    variancesInitializer->setTomoInputs(m_TomoInputs);
    variancesInitializer->setGeometry(m_Geometry);
    variancesInitializer->setObservers(getObservers());
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
    //  std::cout << "------------Initial Variance-----------" << std::endl;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_Sinogram->InitialVariance->d[i_theta] = 1;
      // std::cout << "Tilt: " << i_theta << "  Variance: " << sinogram->InitialVariance->d[i_theta] << std::endl;
    }
  }

  return 0;
}


// -----------------------------------------------------------------------------
// Initialize the Geometry data from a rough reconstruction
// -----------------------------------------------------------------------------
int SOCEngine::initializeRoughReconstructionData()
{
  InitialReconstructionInitializer::Pointer geomInitializer = InitialReconstructionInitializer::NullPointer();
  std::string extension = MXAFileInfo::extension(m_TomoInputs->initialReconFile);
  if (m_TomoInputs->initialReconFile.empty() == true)
  {
    // This will just initialize all the values to Zero (0)
    geomInitializer = InitialReconstructionInitializer::New();
  }
  else if (extension.compare("bin") == 0 )
  {
    // This will read the values from a binary file
    geomInitializer = InitialReconstructionBinReader::NewInitialReconstructionInitializer();
  }
  else if (extension.compare(".mrc") == 0)
  {
    notify("We are not dealing with mrc volume files. The program will now end.", 0, Observable::UpdateErrorMessage);
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
  geomInitializer->setObservers(getObservers());
  geomInitializer->execute();

  if(geomInitializer->getErrorCondition() < 0)
  {
    notify("Error reading Initial Reconstruction Data from File", 100, Observable::UpdateProgressValueAndMessage);
    setErrorCondition(geomInitializer->getErrorCondition());
    return -1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::initializeROIMask(UInt8ImageType::Pointer Mask)
{
  DATA_TYPE x = 0.0;
  DATA_TYPE z = 0.0;
  for (uint16_t i = 0; i < m_Geometry->N_z; i++)
  {
    for (uint16_t j = 0; j < m_Geometry->N_x; j++)
    {
      x = m_Geometry->x0 + ((DATA_TYPE)j + 0.5) * m_TomoInputs->delta_xz;
      z = m_Geometry->z0 + ((DATA_TYPE)i + 0.5) * m_TomoInputs->delta_xz;
      if(x >= -(m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && x <= (m_Sinogram->N_r * m_Sinogram->delta_r) / 2 && z >= -m_TomoInputs->LengthZ / 2
          && z <= m_TomoInputs->LengthZ / 2)
      {
        Mask->d[i][j] = 1;
      }
      else
      {
        Mask->d[i][j] = 0;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::gainAndOffsetInitialization(ScaleOffsetParamsPtr NuisanceParams)
{
  DATA_TYPE sum = 0;
  DATA_TYPE temp = 0;
  for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
  {
    // Gains
    NuisanceParams->I_0->d[k] = m_Sinogram->InitialGain->d[k];
    // Offsets
    NuisanceParams->mu->d[k] = m_Sinogram->InitialOffset->d[k];
#ifdef GEOMETRIC_MEAN_CONSTRAINT
    sum += log(NuisanceParams->I_0->d[k]);
#else
    sum += NuisanceParams->I_0->d[k];
#endif
  }
  sum /= m_Sinogram->N_theta;

#ifdef GEOMETRIC_MEAN_CONSTRAINT
  sum=exp(sum);
  printf("The geometric mean of the gains is %lf\n",sum);

  //Checking if the input parameters satisfy the target geometric mean
  if(fabs(sum - m_Sinogram->targetGain) > 1e-5)
  {
    printf("The input paramters dont meet the constraint..Renormalizing\n");
    temp = exp(log(m_Sinogram->targetGain) - log(sum));
    for (k = 0; k < m_Sinogram->N_theta; k++)
    {
      NuisanceParams->I_0->d[k]=m_Sinogram->InitialGain->d[k]*temp;
    }
  }
#else
  printf("The Arithmetic mean of the constraint is %lf\n", sum);
  if(sum - m_Sinogram->targetGain > 1e-5)
  {
    printf("Arithmetic Mean Constraint not met..renormalizing\n");
    temp = m_Sinogram->targetGain / sum;
    for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
    {
      NuisanceParams->I_0->d[k] = m_Sinogram->InitialGain->d[k] * temp;
    }
  }
#endif

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::initializeHt(RealVolumeType::Pointer H_t)
{
  DATA_TYPE ProfileCenterT;
  for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
  {
    for (int i = 0; i < DETECTOR_RESPONSE_BINS; i++)
    {
      ProfileCenterT = i * OffsetT;
      if(m_TomoInputs->delta_xy >= m_Sinogram->delta_t)
      {
        if(ProfileCenterT <= ((m_TomoInputs->delta_xy / 2) - (m_Sinogram->delta_t / 2)))
        {
          H_t->d[0][k][i] = m_Sinogram->delta_t;
        }
        else
        {
          H_t->d[0][k][i] = -1 * ProfileCenterT + (m_TomoInputs->delta_xy / 2) + m_Sinogram->delta_t / 2;
        }
        if(H_t->d[0][k][i] < 0)
        {
          H_t->d[0][k][i] = 0;
        }

      }
      else
      {
        if(ProfileCenterT <= m_Sinogram->delta_t / 2 - m_TomoInputs->delta_xy / 2)
        {
          H_t->d[0][k][i] = m_TomoInputs->delta_xy;
        }
        else
        {
          H_t->d[0][k][i] = -ProfileCenterT + (m_TomoInputs->delta_xy / 2) + m_Sinogram->delta_t / 2;
        }

        if(H_t->d[0][k][i] < 0)
        {
          H_t->d[0][k][i] = 0;
        }

      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::initializeVolume(RealVolumeType::Pointer Y_Est, double value)
{
  for (uint16_t i = 0; i < m_Sinogram->N_theta; i++)
  {
    for (uint16_t j = 0; j < m_Sinogram->N_r; j++)
    {
      for (uint16_t k = 0; k < m_Sinogram->N_t; k++)
      {
        Y_Est->d[i][j][k] = value;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::storeVoxelResponse(RealVolumeType::Pointer H_t,  AMatrixCol* VoxelLineResponse)
{
  DATA_TYPE ProfileThickness = 0.0;
  DATA_TYPE y = 0.0;
  DATA_TYPE t = 0.0;
  DATA_TYPE tmin;
  DATA_TYPE tmax;
  int16_t slice_index_min, slice_index_max;
  DATA_TYPE center_t,delta_t;
  int16_t index_delta_t;
  DATA_TYPE w3,w4;

  //Storing the response along t-direction for each voxel line
  notify("Storing the response along Y-direction for each voxel line", 0, Observable::UpdateProgressMessage);
  for (uint16_t i =0; i < m_Geometry->N_y; i++)
  {
    y = ((DATA_TYPE)i + 0.5) * m_TomoInputs->delta_xy + m_Geometry->y0;
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
      center_t = ((DATA_TYPE)i_t + 0.5) * m_Sinogram->delta_t + m_Sinogram->T0;
      delta_t = fabs(center_t - t);
      index_delta_t = static_cast<uint16_t>(floor(delta_t / OffsetT));
      if(index_delta_t < DETECTOR_RESPONSE_BINS)
      {
        w3 = delta_t - (DATA_TYPE)(index_delta_t) * OffsetT;
        w4 = ((DATA_TYPE)index_delta_t + 1) * OffsetT - delta_t;
        ProfileThickness = (w4 / OffsetT) * H_t->d[0][0][index_delta_t]
            + (w3 / OffsetT) * H_t->d[0][0][index_delta_t + 1 < DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : DETECTOR_RESPONSE_BINS - 1];
    //  ProfileThickness = (w4 / OffsetT) * detectorResponse->d[0][uint16_t(floor(m_Sinogram->N_theta/2))][index_delta_t]
    //  + (w3 / OffsetT) * detectorResponse->d[0][uint16_t(floor(m_Sinogram->N_theta/2))][index_delta_t + 1 < DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : DETECTOR_RESPONSE_BINS - 1];
    }
      else
      {
        ProfileThickness = 0;
      }

      if(ProfileThickness != 0) //Store the response of this slice
      {
#ifdef DEBUG
        printf("%d %lf\n", i_t, ProfileThickness);
#endif
        VoxelLineResponse[i].values[VoxelLineResponse[i].count] = ProfileThickness;
        VoxelLineResponse[i].index[VoxelLineResponse[i].count++] = i_t;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Calculate Error Sinogram
// Also compute weights of the diagonal covariance matrix
// -----------------------------------------------------------------------------
void SOCEngine::calculateMeasurementWeight(RealVolumeType::Pointer Weight,
                                           ScaleOffsetParamsPtr NuisanceParams,
                                           RealVolumeType::Pointer ErrorSino,
                                           RealVolumeType::Pointer Y_Est)
{
  DATA_TYPE checksum = 0;
  START_TIMER;
  for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++) //slice index
  {
#ifdef NOISE_MODEL
    {
      NuisanceParams->alpha->d[i_theta] = m_Sinogram->InitialVariance->d[i_theta]; //Initialize the refinement parameters from any previous run
    }
#endif
    checksum = 0;
    for (int16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        ErrorSino->d[i_theta][i_r][i_t] = m_Sinogram->counts->d[i_theta][i_r][i_t] - Y_Est->d[i_theta][i_r][i_t] - NuisanceParams->mu->d[i_theta];

        if(m_Sinogram->counts->d[i_theta][i_r][i_t] != 0)
        {
          Weight->d[i_theta][i_r][i_t] = 1.0 / m_Sinogram->counts->d[i_theta][i_r][i_t];
#ifdef NOISE_MODEL
          {
            Weight->d[i_theta][i_r][i_t] /= NuisanceParams->alpha->d[i_theta];
          }
#endif
        }
        else
        {
          Weight->d[i_theta][i_r][i_t] = 0;
        }

#ifdef FORWARD_PROJECT_MODE
        temp=Y_Est->d[i_theta][i_r][i_t]/NuisanceParams->I_0->d[i_theta];
        fwrite(&temp,sizeof(DATA_TYPE),1,Fp6);
#endif
        if(Weight->d[i_theta][i_r][i_t] < 0)
        {
          std::cout << m_Sinogram->counts->d[i_theta][i_r][i_t] << "    " << NuisanceParams->alpha->d[i_theta] << std::endl;
        }

        checksum += Weight->d[i_theta][i_r][i_t];
      }
    }
#ifdef DEBUG
    printf("Check sum of Diagonal Covariance Matrix= %lf\n", checksum);
#endif
  }
  STOP_TIMER;
  std::cout << std::endl;
  std::string indent(" ");
  PRINT_TIME("Computing Weights");
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int SOCEngine::calculateCost(CostData::Pointer cost,
                             RealVolumeType::Pointer Weight,
                             RealVolumeType::Pointer ErrorSino)
{
  DATA_TYPE cost_value = computeCost(ErrorSino, Weight);
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
void SOCEngine::updateWeights(RealVolumeType::Pointer Weight,
                              ScaleOffsetParamsPtr NuisanceParams,
                              RealVolumeType::Pointer ErrorSino)
{
  DATA_TYPE AverageVarUpdate = 0; //absolute sum of the gain updates
  DATA_TYPE AverageMagVar = 0; //absolute sum of the initial gains
  DATA_TYPE sum = 0;

  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
  {
    uint32_t NumNonZeroEntries = 0;
    sum = 0;
    //Factoring out the variance parameter from the Weight matrix
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        if(m_Sinogram->counts->d[i_theta][i_r][i_t] != 0)
        {
          Weight->d[i_theta][i_r][i_t] = 1.0 / m_Sinogram->counts->d[i_theta][i_r][i_t];
          NumNonZeroEntries++;
        }
      }
    }

    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        sum += (ErrorSino->d[i_theta][i_r][i_t] * ErrorSino->d[i_theta][i_r][i_t] * Weight->d[i_theta][i_r][i_t]); //Changed to only account for the counts
      }
    }
    sum /= NumNonZeroEntries; //(m_Sinogram->N_r*m_Sinogram->N_t);

    AverageMagVar += fabs(NuisanceParams->alpha->d[i_theta]);
    AverageVarUpdate += fabs(sum - NuisanceParams->alpha->d[i_theta]);
    NuisanceParams->alpha->d[i_theta] = sum;
    //Update the weight for ICD updates
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
    {
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        if(NuisanceParams->alpha->d[i_theta] != 0 && m_Sinogram->counts->d[i_theta][i_r][i_t] != 0)
        {
          Weight->d[i_theta][i_r][i_t] = 1.0 / (m_Sinogram->counts->d[i_theta][i_r][i_t] * NuisanceParams->alpha->d[i_theta]);
        }
        else
        {
          Weight->d[i_theta][i_r][i_t] = 0;
        }

      }
    }

  }

#ifdef DEBUG
  std::cout << "Noise Model Weights:" << std::endl;
  std::cout << "Tilt\tWeight" << std::endl;
  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
  {
    std::cout << i_theta << "\t" << NuisanceParams->alpha->d[i_theta] << std::endl;
  }
#endif
  DATA_TYPE VarRatio = AverageVarUpdate / AverageMagVar;
  std::cout << "Ratio of change in Variance " << VarRatio << std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::writeNuisanceParameters(ScaleOffsetParamsPtr NuisanceParams)
{
  NuisanceParamWriter::Pointer nuisanceBinWriter = NuisanceParamWriter::New();
  nuisanceBinWriter->setSinogram(m_Sinogram);
  nuisanceBinWriter->setTomoInputs(m_TomoInputs);
  nuisanceBinWriter->setObservers(getObservers());
  nuisanceBinWriter->setNuisanceParams(NuisanceParams.get());
#ifdef JOINT_ESTIMATION
  nuisanceBinWriter->setFileName(m_TomoInputs->gainsOutputFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_I_O);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
  }

  nuisanceBinWriter->setFileName(m_TomoInputs->offsetsOutputFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_mu);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
   setErrorCondition(-1);
   notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
  }
#endif

#ifdef NOISE_MODEL
  nuisanceBinWriter->setFileName(m_TomoInputs->varianceOutputFile);
  nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_alpha);
  nuisanceBinWriter->execute();
  if (nuisanceBinWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
  }
#endif//Noise Model

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::writeSinogramFile(ScaleOffsetParamsPtr NuisanceParams, RealVolumeType::Pointer Final_Sinogram)
{
  // Write the Sinogram out to a file
  SinogramBinWriter::Pointer sinogramWriter = SinogramBinWriter::New();
  sinogramWriter->setSinogram(m_Sinogram);
  sinogramWriter->setTomoInputs(m_TomoInputs);
  sinogramWriter->setObservers(getObservers());
  sinogramWriter->setNuisanceParams(NuisanceParams);
  sinogramWriter->setData(Final_Sinogram);
  sinogramWriter->execute();
  if (sinogramWriter->getErrorCondition() < 0)
  {
    setErrorCondition(-1);
    notify(sinogramWriter->getErrorMessage().c_str(), 100, Observable::UpdateProgressValueAndMessage);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::writeReconstructionFile()
{
  // Write the Reconstruction out to a file
  RawGeometryWriter::Pointer writer = RawGeometryWriter::New();
  writer->setGeometry(m_Geometry);
  writer->setFilePath(m_TomoInputs->reconstructedOutputFile);
  writer->setObservers(getObservers());
  writer->execute();
  if (writer->getErrorCondition() < 0)
  {
    setErrorCondition(writer->getErrorCondition());
    notify("Error Writing the Raw Geometry", 100, Observable::UpdateProgressValueAndMessage);
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::writeVtkFile()
{
  std::string vtkFile(m_TomoInputs->tempDir);
  vtkFile = vtkFile.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::ReconstructedVtkFile);

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

  VtkScalarWriter* w0 = static_cast<VtkScalarWriter*>(new TomoOutputScalarWriter(m_Geometry.get()));
  w0->setWriteBinaryFiles(true);
  scalarsToWrite.push_back(w0);

  int error = vtkWriter.write<DimsAndRes>(vtkFile, &dimsAndRes, scalarsToWrite);
  if (error < 0)
  {
    std::cout << "Error writing vtk file '" << vtkFile << "'" << std::endl;
  }
  delete w0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SOCEngine::writeMRCFile()
{
  /* Write the output to the MRC File */
   std::string mrcFile (m_TomoInputs->tempDir);
   mrcFile = mrcFile.append(MXADir::getSeparator()).append(ScaleOffsetCorrection::ReconstructedMrcFile);
   MRCWriter::Pointer mrcWriter = MRCWriter::New();
   mrcWriter->setOutputFile(mrcFile);
   mrcWriter->setGeometry(m_Geometry);
   mrcWriter->write();
}
