
#include "HAADF_ForwardModel.h"

#include "MBIRLib/IOFilters/MRCWriter.h"
#include "MBIRLib/IOFilters/RawGeometryWriter.h"
#include "MBIRLib/IOFilters/NuisanceParamWriter.h"
#include "MBIRLib/IOFilters/SinogramBinWriter.h"
#include "MBIRLib/IOFilters/VTKFileWriters.hpp"
#include "MBIRLib/IOFilters/AvizoUniformCoordinateWriter.h"

#include "MBIRLib/IOFilters/NuisanceParamReader.h"
#include "MBIRLib/IOFilters/GainsOffsetsReader.h"

#include "MBIRLib/GenericFilters/ComputeInitialOffsets.h"



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADF_ForwardModel::HAADF_ForwardModel()
{}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
HAADF_ForwardModel::~HAADF_ForwardModel()
{}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ForwardModel::writeNuisanceParameters(SinogramPtr sinogram)
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
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateErrorMessage);
    }

    nuisanceBinWriter->setFileName(m_TomoInputs->offsetsOutputFile);
    nuisanceBinWriter->setDataToWrite(NuisanceParamWriter::Nuisance_mu);
    nuisanceBinWriter->setData(m_Mu);
    nuisanceBinWriter->execute();
    if(nuisanceBinWriter->getErrorCondition() < 0)
    {
      setErrorCondition(-1);
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateErrorMessage);
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
      notify(nuisanceBinWriter->getErrorMessage().c_str(), 100, Observable::UpdateErrorMessage);
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
void HAADF_ForwardModel::writeSinogramFile(SinogramPtr sinogram, RealVolumeType::Pointer finalSinogram)
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
    notify(sinogramWriter->getErrorMessage().c_str(), 100, Observable::UpdateErrorMessage);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ForwardModel::writeReconstructionFile(const std::string& filepath)
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
    notify("Error Writing the Raw Geometry", 100, Observable::UpdateErrorMessage);
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ForwardModel::writeVtkFile(const std::string& vtkFile, uint16_t cropStart, uint16_t cropEnd)
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
void HAADF_ForwardModel::writeMRCFile(const std::string& mrcFile, uint16_t cropStart, uint16_t cropEnd)
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
void HAADF_ForwardModel::writeAvizoFile(const std::string& file, uint16_t cropStart, uint16_t cropEnd)
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


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ForwardModel::createInitialGainsData()
{
  std::stringstream ss;
  /* ********************* Initialize the Gains Array **************************/
  size_t gains_dims[1] =
  { m_Sinogram->N_theta };
  m_InitialGain = RealArrayType::New(gains_dims, "sinogram->InitialGain");
  if(m_TomoInputs->gainsInputFile.empty() == false)
  {
    // Read the initial Gains from a File
    NuisanceParamReader::Pointer gainsInitializer = NuisanceParamReader::New();
    gainsInitializer->setFileName(m_TomoInputs->gainsInputFile);
    gainsInitializer->setData(m_InitialGain);
    gainsInitializer->setSinogram(m_Sinogram);
    gainsInitializer->setAdvParams(m_AdvParams);
    gainsInitializer->setTomoInputs(m_TomoInputs);
    gainsInitializer->setGeometry(m_Geometry);
    gainsInitializer->setObservers(getObservers());
    gainsInitializer->execute();
    if(gainsInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Gains from Data file", 100, Observable::UpdateErrorMessage);
      setErrorCondition(gainsInitializer->getErrorCondition());
      return -1;
    }
  }
  else
  {
    // Set the values to the target gain value set by the user
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_InitialGain->d[i_theta] = m_TomoInputs->targetGain;
    }
  }
  /********************REMOVE************************/
  ss << "HARD WIRED TARGET GAIN" << std::endl;
  m_TargetGain = m_TomoInputs->targetGain; //TARGET_GAIN;
  ss << "Target Gain: " << m_TargetGain << std::endl;
  /*************************************************/

  if (getVeryVerbose())
  {
    std::cout << ss.str() << std::endl;
  }
  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ForwardModel::createInitialOffsetsData()
{
  std::stringstream ss;
  /* ********************* Initialize the Offsets Array **************************/
  size_t offsets_dims[1] =
  { m_Sinogram->N_theta };
  m_InitialOffset = RealArrayType::New(offsets_dims, "sinogram->InitialOffset");
  if(m_TomoInputs->offsetsInputFile.empty() == false)
  {
    // Read the initial offsets from a File
    NuisanceParamReader::Pointer offsetsInitializer = NuisanceParamReader::New();
    offsetsInitializer->setFileName(m_TomoInputs->offsetsInputFile);
    offsetsInitializer->setData(m_InitialOffset);
    offsetsInitializer->setSinogram(m_Sinogram);
    offsetsInitializer->setAdvParams(m_AdvParams);
    offsetsInitializer->setTomoInputs(m_TomoInputs);
    offsetsInitializer->setGeometry(m_Geometry);
    offsetsInitializer->setObservers(getObservers());
    offsetsInitializer->execute();
    if(offsetsInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Offsets from Data file", 100, Observable::UpdateErrorMessage);
      setErrorCondition(offsetsInitializer->getErrorCondition());
      return -1;
    }
  }
  else if(m_TomoInputs->useDefaultOffset == true)
  {
    // Set the values to the default offset value set by the user
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_InitialOffset->d[i_theta] = m_TomoInputs->defaultOffset;
    }
  }
  else
  {
    // Compute the initial offset values from the data
    ComputeInitialOffsets::Pointer initializer = ComputeInitialOffsets::New();
    initializer->setSinogram(m_Sinogram);
    initializer->setTomoInputs(m_TomoInputs);
    initializer->setAdvParams(m_AdvParams);
    initializer->setInitialOffset(m_InitialOffset);
    initializer->setObservers(getObservers());
    initializer->setVerbose(getVerbose());
    initializer->setVeryVerbose(getVeryVerbose());
    initializer->execute();
    if(initializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Offsets Data file", 100, Observable::UpdateErrorMessage);
      setErrorCondition(initializer->getErrorCondition());
      return -1;
    }
  }
  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int HAADF_ForwardModel::createInitialVariancesData()
{

  /* ********************* Initialize the Variances Array **************************/
  size_t variance_dims[1] =
  { m_Sinogram->N_theta };
  m_InitialVariance = RealArrayType::New(variance_dims, "sinogram->InitialVariance");
  if(m_TomoInputs->varianceInputFile.empty() == false)
  {
    // Read the initial variances from a File
    NuisanceParamReader::Pointer variancesInitializer = NuisanceParamReader::New();
    variancesInitializer->setFileName(m_TomoInputs->varianceInputFile);
    variancesInitializer->setData(m_InitialVariance);
    variancesInitializer->setSinogram(m_Sinogram);
    variancesInitializer->setTomoInputs(m_TomoInputs);
    variancesInitializer->setGeometry(m_Geometry);
    variancesInitializer->setAdvParams(m_AdvParams);
    variancesInitializer->setObservers(getObservers());
    variancesInitializer->setVerbose(getVerbose());
    variancesInitializer->setVeryVerbose(getVeryVerbose());
    variancesInitializer->execute();
    if(variancesInitializer->getErrorCondition() < 0)
    {
      notify("Error initializing Input Variances from Data file", 100, Observable::UpdateErrorMessage);
      setErrorCondition(variancesInitializer->getErrorCondition());
      return -1;
    }
  }
  else
  {
    std::stringstream ss;
    for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    {
      m_InitialVariance->d[i_theta] = m_TomoInputs->defaultVariance;
      ss << "Tilt: " << i_theta << "  Variance: " << m_InitialVariance->d[i_theta] << std::endl;
    }
    if (getVeryVerbose())
    {
      std::cout << ss.str() << std::endl;
    }
  }

  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ForwardModel::initializeHt(RealVolumeType::Pointer H_t, Real_t OffsetT)
{
  Real_t ProfileCenterT;
  for (uint16_t k = 0; k < m_Sinogram->N_theta; k++)
  {
    for (unsigned int i = 0; i < m_AdvParams->DETECTOR_RESPONSE_BINS; i++)
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
void HAADF_ForwardModel::selectorInitialization(size_t dims[3])
{
  //This variable selects which entries to retain in the sinogram
  m_BraggSelector = UInt8VolumeType::New(dims, "Selector");

}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void HAADF_ForwardModel::computeBraggSelector(RealVolumeType::Pointer ErrorSino, RealVolumeType::Pointer Weight)
{
  Real_t sum = 0;
  for (uint16_t i_theta = 0; i_theta < m_Sinogram->N_theta; i_theta++)
    for (uint16_t i_r = 0; i_r < m_Sinogram->N_r; i_r++)
      for (uint16_t i_t = 0; i_t < m_Sinogram->N_t; i_t++)
      {
        size_t idx = Weight->calcIndex(i_theta, i_r, i_t);
        if(ErrorSino->d[idx] * ErrorSino->d[idx] * Weight->d[idx] < m_BraggThreshold * m_BraggThreshold)
	  { m_BraggSelector->d[idx] = 1; sum+=1;}
        else
        { m_BraggSelector->d[idx] = 0; }
      }
  std::cout<<"Non-zero selector entries = "<<sum<<std::endl;
}

