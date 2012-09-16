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

#include "ReconstructionEngine.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/GenericFilters/MRCSinogramInitializer.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Common/EIMTime.h"
#include "MBIRLib/HAADF/UpdateYSlice.h"


#define START_TIMER uint64_t startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER uint64_t stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME(msg)\
std::cout << indent << msg << ": " << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

// Read the Input data from the supplied data file
// We are scoping here so the various readers are automatically cleaned up before
// the code goes any farther
int ReconstructionEngine::readInputData()
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
  dataReader->setAdvParams(m_AdvParams);
  dataReader->setObservers(getObservers());
  dataReader->setVerbose(getVerbose());
  dataReader->setVeryVerbose(getVeryVerbose());
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
// Initialize the Geometry data from a rough reconstruction
// -----------------------------------------------------------------------------
int ReconstructionEngine::initializeRoughReconstructionData()
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
    notify("Error reading Initial Reconstruction Data from File", 100, Observable::UpdateProgressValueAndMessage);
    setErrorCondition(geomInitializer->getErrorCondition());
    return -1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::computeOriginalXDims(uint16_t &cropStart, uint16_t &cropEnd)
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
void ReconstructionEngine::initializeHt(RealVolumeType::Pointer H_t, Real_t OffsetT)
{
  Real_t ProfileCenterT;
  for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
  {
    for (unsigned int i = 0; i <m_AdvParams->DETECTOR_RESPONSE_BINS; i++)
    {
      ProfileCenterT = i * OffsetT;
      if(m_TomoInputs->delta_xy >= m_Sinogram->delta_t)
      {
        if(ProfileCenterT <= ((m_TomoInputs->delta_xy / 2) - (m_Sinogram->delta_t / 2)))
        {
          H_t->setValue(m_Sinogram->delta_t, 0, k, i);
        }
        else
        {
          H_t->setValue(-1 * ProfileCenterT + (m_TomoInputs->delta_xy / 2) + m_Sinogram->delta_t / 2, 0, k, i);
        }
        if(H_t->getValue(0, k, i) < 0)
        {
          H_t->setValue(0, 0, k, i);
        }

      }
      else
      {
        if(ProfileCenterT <= m_Sinogram->delta_t / 2 - m_TomoInputs->delta_xy / 2)
        {
          H_t->setValue(m_TomoInputs->delta_xy, 0, k, i);
        }
        else
        {
          H_t->setValue(-ProfileCenterT + (m_TomoInputs->delta_xy / 2) + m_Sinogram->delta_t / 2, 0, k, i);
        }

        if(H_t->getValue(0, k, i) < 0)
        {
          H_t->setValue(0, 0, k, i);
        }

      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::initializeVolume(RealVolumeType::Pointer Y_Est, double value)
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
void ReconstructionEngine::storeVoxelResponse(RealVolumeType::Pointer H_t,
                                              std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
                                              HAADFDetectorParameters::Pointer haadfParameters)
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

  const Real_t offsetT = haadfParameters->getOffsetT();
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
      index_delta_t = static_cast<uint16_t>(floor(delta_t / offsetT));
      if(index_delta_t < m_AdvParams->DETECTOR_RESPONSE_BINS)
      {
        w3 = delta_t - (Real_t)(index_delta_t) * offsetT;
        w4 = ((Real_t)index_delta_t + 1) * offsetT - delta_t;
        uint16_t ttmp = index_delta_t + 1 < m_AdvParams->DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : m_AdvParams->DETECTOR_RESPONSE_BINS - 1;
        ProfileThickness = (w4 / offsetT) * H_t->getValue(0, 0, index_delta_t) + (w3 / offsetT) * H_t->getValue(0, 0, ttmp);
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
        HAADFAMatrixCol::Pointer vlr = VoxelLineResponse[i];
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
void ReconstructionEngine::calculateArithmeticMean()
{

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::writeReconstructionFile(const std::string &filepath)
{
  // Write the Reconstruction out to a file
  RawGeometryWriter::Pointer writer = RawGeometryWriter::New();
  writer->setGeometry(m_Geometry);
  writer->setFilePath(filepath);
  writer->setAdvParams(m_AdvParams);
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
void ReconstructionEngine::writeVtkFile(const std::string &vtkFile, uint16_t cropStart, uint16_t cropEnd)
{
  std::stringstream ss;
  ss << "Writing VTK file to '" << vtkFile << "'";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);

  VTKStructuredPointsFileWriter vtkWriter;
  vtkWriter.setWriteBinaryFiles(true);


  DimsAndRes dimsAndRes;
  dimsAndRes.xStart = cropStart;
  dimsAndRes.xEnd = cropEnd;
  dimsAndRes.yStart = 0;
  dimsAndRes.yEnd = m_Geometry->N_y;
  dimsAndRes.zStart = 0;
  dimsAndRes.zEnd = m_Geometry->N_z;
  dimsAndRes.resx = 1.0f;
  dimsAndRes.resy = 1.0f;
  dimsAndRes.resz = 1.0f;

  std::vector<VtkScalarWriter*> scalarsToWrite;

  VtkScalarWriter* w0 = static_cast<VtkScalarWriter*>(new TomoOutputScalarWriter(m_Geometry.get()));
  w0->setXDims(cropStart, cropEnd);
  w0->setYDims(0, m_Geometry->N_y);
  w0->setZDims(0, m_Geometry->N_z);
  w0->setWriteBinaryFiles(true);
  scalarsToWrite.push_back(w0);

  int error = vtkWriter.write<DimsAndRes>(vtkFile, &dimsAndRes, scalarsToWrite);
  if(error < 0)
  {
    ss.str("");
    ss << "Error writing vtk file\n    '" << vtkFile << "'" << std::endl;
    setErrorCondition(-12);
    notify(ss.str(), 0, Observable::UpdateErrorMessage);
  }
  delete w0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::writeMRCFile(const std::string &mrcFile, uint16_t cropStart, uint16_t cropEnd)
{
  /* Write the output to the MRC File */
  std::stringstream ss;
  ss.str("");
  ss << "Writing MRC file to '" << mrcFile << "'";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);

  MRCWriter::Pointer mrcWriter = MRCWriter::New();
  mrcWriter->setOutputFile(mrcFile);
  mrcWriter->setGeometry(m_Geometry);
  mrcWriter->setAdvParams(m_AdvParams);
  mrcWriter->setXDims(cropStart, cropEnd);
  mrcWriter->setYDims(0, m_Geometry->N_y);
  mrcWriter->setZDims(0, m_Geometry->N_z);
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

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::writeAvizoFile(const std::string &file, uint16_t cropStart, uint16_t cropEnd)
{
  //  std::cout << "Writing Avizo file " << file << std::endl;

  /* Write the output to the Avizo File */
  std::stringstream ss;
  ss.str("");
  ss << "Writing Avizo file to '" << file << "'";
  notify(ss.str(), 0, Observable::UpdateProgressMessage);

  AvizoUniformCoordinateWriter::Pointer writer = AvizoUniformCoordinateWriter::New();
  writer->setOutputFile(file);
  writer->setGeometry(m_Geometry);
  writer->setAdvParams(m_AdvParams);
  writer->setTomoInputs(m_TomoInputs);
  writer->setXDims(cropStart, cropEnd);
  writer->setYDims(0, m_Geometry->N_y);
  writer->setZDims(0, m_Geometry->N_z);
  writer->setObservers(getObservers());
  writer->setWriteBinaryFile(true);
  writer->execute();
  if(writer->getErrorCondition() < 0)
  {
    ss.str("");
    ss << "Error writing Avizo file\n    '" << file << "'" << std::endl;
    setErrorCondition(writer->getErrorCondition());
    notify(ss.str(), 0, Observable::UpdateErrorMessage);
  }

}


uint8_t ReconstructionEngine::updateVoxels( //SinogramPtr sinogram,
                                           //GeometryPtr geometry,
                                           int16_t OuterIter,
                                           int16_t Iter,
                                           UInt8Image_t::Pointer VisitCount,
                                           std::vector<HAADFAMatrixCol::Pointer> &TempCol,
                                           RealVolumeType::Pointer ErrorSino,
                                           std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
                                           CostData::Pointer cost,
                                           QGGMRF::QGGMRF_Values* qggmrf_values)
{
    size_t dims[3];
    dims[0] = m_Geometry->N_z; //height
    dims[1] = m_Geometry->N_x; //width
    dims[2] = 0;

    RealImageType::Pointer magUpdateMap = RealImageType::New(dims, "Update Map for voxel lines");
    RealImageType::Pointer filtMagUpdateMap = RealImageType::New(dims, "Filter Update Map for voxel lines");
    UInt8Image_t::Pointer magUpdateMask = UInt8Image_t::New(dims, "Update Mask for selecting voxel lines NHICD");

#if ROI
    UInt8Image_t::Pointer mask;
    dims[0] = m_Geometry->N_z;
    dims[1] = m_Geometry->N_x;
    mask = UInt8Image_t::New(dims, "Mask");
    initializeROIMask(mask);
#endif

    unsigned int updateType = VoxelUpdateType::RegularRandomOrderUpdate;
#ifdef NHICD
    if(0 == reconInnerIter % 2)
    {
        updateType = VoxelUpdateType::HomogeniousUpdate;
    }
    else
    {
        updateType = VoxelUpdateType::NonHomogeniousUpdate;
    }
#else

#endif//NHICD end if

#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
    tbb::task_scheduler_init init;
    int m_NumThreads = init.default_num_threads();
#else
    int m_NumThreads = 1;
#endif

    std::stringstream ss;
    uint8_t exit_status = 1; //Indicates normal exit ; else indicates to stop inner iterations
    uint16_t subIterations = 1;
    std::string indent("    ");
    uint8_t err = 0;

    if(updateType == VoxelUpdateType::RegularRandomOrderUpdate)
    {
        ss << indent << "Regular Random Order update of Voxels" << std::endl;
    }
    else if(updateType == VoxelUpdateType::HomogeniousUpdate)
    {
        ss << indent << "Homogenous update of voxels" << std::endl;
    }
    else if(updateType == VoxelUpdateType::NonHomogeniousUpdate)
    {
        ss << indent << "Non Homogenous update of voxels" << std::endl;
        subIterations = NUM_NON_HOMOGENOUS_ITER;
    }
    else
    {
        ss << indent << "Unknown Voxel Update Type. Returning Now" << std::endl;
        notify(ss.str(), 0, Observable::UpdateErrorMessage);
        return exit_status;
    }

    if(getVerbose())
    {
        std::cout << ss.str() << std::endl;
    }

    Real_t NH_Threshold = 0.0;
    int totalLoops = m_TomoInputs->NumOuterIter * m_TomoInputs->NumIter;

    for (uint16_t NH_Iter = 0; NH_Iter < subIterations; ++NH_Iter)
    {
        ss.str("");
        ss << "Outer Iteration: " << OuterIter << " of " << m_TomoInputs->NumOuterIter;
        ss << "   Inner Iteration: " << Iter << " of " << m_TomoInputs->NumIter;
        ss << "   SubLoop: " << NH_Iter << " of " << subIterations;
        float currentLoop = static_cast<float>(OuterIter * m_TomoInputs->NumIter + Iter);
        notify(ss.str(), currentLoop / totalLoops * 100.0f, Observable::UpdateProgressValueAndMessage);
        if(updateType == VoxelUpdateType::NonHomogeniousUpdate)
        {
            //Compute VSC and create a map of pixels that are above the threshold value
            ComputeVSC(magUpdateMap, filtMagUpdateMap);
            START_TIMER;
            NH_Threshold = SetNonHomThreshold(magUpdateMap);
            STOP_TIMER;
            PRINT_TIME("  SetNonHomThreshold");
            std::cout << indent << "NHICD Threshold: " << NH_Threshold << std::endl;
            //Use  filtMagUpdateMap  to find MagnitudeUpdateMask
            //std::cout << "Completed Calculation of filtered magnitude" << std::endl;
            //Calculate the threshold for the top ? % of voxel updates
        }

        //printf("Iter %d\n",Iter);
#if ROI
        //variables used to stop the process
        Real_t AverageUpdate = 0;
        Real_t AverageMagnitudeOfRecon = 0;
#endif

        START_TIMER;
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
        std::vector<int> yCount(m_NumThreads, 0);
        int t = 0;
        for (int y = 0; y < m_Geometry->N_y; ++y)
        {
            yCount[t]++;
            ++t;
            if(t == m_NumThreads)
            {
                t = 0;
            }
        }

        uint16_t yStart = 0;
        uint16_t yStop = 0;

        tbb::task_list taskList;
        Real_t* averageUpdate = (Real_t*)(malloc(sizeof(Real_t) * m_NumThreads));
        ::memset(averageUpdate, 0, sizeof(Real_t) * m_NumThreads);
        Real_t* averageMagnitudeOfRecon = (Real_t*)(malloc(sizeof(Real_t) * m_NumThreads));
        ::memset(averageMagnitudeOfRecon, 0, sizeof(Real_t) * m_NumThreads);
        for (int t = 0; t < m_NumThreads; ++t)
        {
            yStart = yStop;
            yStop = yStart + yCount[t];
            if(yStart == yStop)
            {
                continue;
            } // Processor has NO tasks to run because we have less Y's than cores

            // std::cout << "Thread: " << t << " yStart: " << yStart << "  yEnd: " << yStop << std::endl;
            UpdateYSlice& a =
            *new (tbb::task::allocate_root()) UpdateYSlice(yStart, yStop, m_Geometry, OuterIter, Iter,
                                                           m_Sinogram, TempCol, ErrorSino,
                                                           VoxelLineResponse, m_ForwardModel.get(), mask,
                                                           magUpdateMap, magUpdateMask, updateType,
                                                           NH_Threshold,
                                                           averageUpdate + t,
                                                           averageMagnitudeOfRecon + t,
                                                           m_AdvParams->ZERO_SKIPPING,
                                                           qggmrf_values);
            taskList.push_back(a);
        }

        tbb::task::spawn_root_and_wait(taskList);
        // Now sum up some values
        for (int t = 0; t < m_NumThreads; ++t)
        {
            AverageUpdate += averageUpdate[t];
            AverageMagnitudeOfRecon += averageMagnitudeOfRecon[t];
        }
        free(averageUpdate);
        free(averageMagnitudeOfRecon);

#else
        uint16_t yStop = geometry->N_y;
        uint16_t yStart = 0;
        UpdateYSlice yVoxelUpdate(yStart, yStop,
                                  geometry,
                                  OuterIter, Iter, sinogram,
                                  m_BFSinogram, TempCol,
                                  ErrorSino, Weight, VoxelLineResponse,
                                  NuisanceParams, Mask,
                                  magUpdateMap, MagUpdateMask,
                                  &m_QGGMRF_Values,
                                  updateType,
                                  NH_Threshold,
                                  &AverageUpdate,
                                  &AverageMagnitudeOfRecon,
                                  m_AdvParams->ZERO_SKIPPING);

        yVoxelUpdate.execute();
#endif
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */STOP_TIMER;
        ss.str("");
        ss << "Inner Iter: " << Iter << " Voxel Update";
        PRINT_TIME(ss.str());


#if ROI
        if(getVerbose())
        {
            std::cout << "Average Update " << AverageUpdate << std::endl;
            std::cout << "Average Mag " << AverageMagnitudeOfRecon << std::endl;
        }
        if(AverageMagnitudeOfRecon > 0)
        {
            if(getVerbose())
            {
                std::cout << Iter + 1 << " " << AverageUpdate / AverageMagnitudeOfRecon << std::endl;
            }
            //Use the stopping criteria if we are performing a full update of all voxels
            if((AverageUpdate / AverageMagnitudeOfRecon) < m_TomoInputs->StopThreshold && updateType != VoxelUpdateType::NonHomogeniousUpdate)
            {
                std::cout << "This is the terminating point " << Iter << std::endl;
                m_TomoInputs->StopThreshold *= m_AdvParams->THRESHOLD_REDUCTION_FACTOR; //Reducing the thresold for subsequent iterations
                std::cout << "New threshold" << m_TomoInputs->StopThreshold << std::endl;
                exit_status = 0;
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
            TempPointer = geometry->Object;
            NumOfBytesWritten=fwrite(&(geometry->Object->d[0][0][0]), sizeof(Real_t),geometry->N_x*geometry->N_y*geometry->N_z, Fp3);
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



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionEngine::initializeROIMask(UInt8Image_t::Pointer Mask)
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
void ReconstructionEngine::ComputeVSC(RealImageType::Pointer magUpdateMap, RealImageType::Pointer filtMagUpdateMap)
{
    Real_t filter_op = 0;
    // int err = 0;
    FILE *Fp = NULL;
    MAKE_OUTPUT_FILE(Fp, m_TomoInputs->tempDir, ScaleOffsetCorrection::MagnitudeMapFile);
    if(errno < 0)
    {

    }
    fwrite(magUpdateMap->getPointer(0, 0), m_Geometry->N_x * m_Geometry->N_z, sizeof(Real_t), Fp);
    fclose(Fp);

    // std::cout<<"Starting to filter the magnitude"<<std::endl;
    // std::cout<<geometry->N_x<<" " <<geometry->N_z<<std::endl;
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
                //		filter_op += k_HammingWindow[p + 2][q + 2] * magUpdateMap->getValue(i + p, j + q);
                    }
                }
            }
            filtMagUpdateMap->setValue(filter_op, i, j);
        }
    }

    for (int16_t i = 0; i < m_Geometry->N_z; i++)
    {
        for (int16_t j = 0; j < m_Geometry->N_x; j++)
        {
            //magUpdateMap->d[i][j]=filtMagUpdateMap->d[i][j];
            magUpdateMap->setValue(filtMagUpdateMap->getValue(i, j), i, j);
        }
    }

    MAKE_OUTPUT_FILE(Fp, m_TomoInputs->tempDir, ScaleOffsetCorrection::FilteredMagMapFile);
    if(errno < 0)
    {

    }
    fwrite(filtMagUpdateMap->getPointer(0, 0), m_Geometry->N_x * m_Geometry->N_z, sizeof(Real_t), Fp);
    fclose(Fp);
}


// -----------------------------------------------------------------------------
// Sort the entries of filtMagUpdateMap and set the threshold to be ? percentile
// -----------------------------------------------------------------------------
Real_t ReconstructionEngine::SetNonHomThreshold(RealImageType::Pointer magUpdateMap)
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
            //TempMagMap->d[i*geometry->N_x+j]=i*geometry->N_x+j;
            TempMagMap->d[i * (uint32_t)m_Geometry->N_x + j] = magUpdateMap->getValue(i, j);
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

    threshold = TempMagMap->d[percentile_index];
    return threshold;
}


