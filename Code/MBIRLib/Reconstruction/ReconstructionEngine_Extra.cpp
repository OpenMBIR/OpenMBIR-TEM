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

#ifdef BF_RECON
void ReconstructionEngine::processRawCounts()
{
    Real_t mean=0;
    for (int16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++) //slice index
    {
        for (int16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
        {
            for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
            {
                size_t counts_idx = m_Sinogram->counts->calcIndex(i_theta, i_r, i_t);
                m_Sinogram->counts->d[counts_idx] += BF_OFFSET;
                m_Sinogram->counts->d[counts_idx] = -log(m_Sinogram->counts->d[counts_idx]/BF_MAX);

                if(m_Sinogram->counts->d[counts_idx] < 0 ) //Clip the log data to be positive
                    m_Sinogram->counts->d[counts_idx] = 0;

                mean+=m_Sinogram->counts->d[counts_idx];
            }
        }
    }
    mean/=(m_Sinogram->N_theta*m_Sinogram->N_r*m_Sinogram->N_t);
    std::cout<<"Mean log value ="<<mean<<std::endl;
}

#endif //BF Recon




