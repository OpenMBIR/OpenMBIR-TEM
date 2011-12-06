/*
 * GainsOffsetsReader.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */
#include "GainsOffsetsReader.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsReader::GainsOffsetsReader() :
m_Inputs(NULL),
m_Sinogram(NULL)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsReader::~GainsOffsetsReader()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsReader::execute()
{
  FILE* Fp = NULL;
  uint16_t view_count = 0;
 double *buffer;
  Fp = fopen(m_Inputs->GainsOffsetsFile.c_str(), "r"); //This file contains the Initial unscatterd counts and background scatter for each view
	//Fp=fopen("/Users/singanallurvenkatakrishnan/Desktop/Work/Tomography/TomoSoftware/HAADFSTEM/C-Code/Data/ConvergedGainOffsetParamsOuter45_AMConstraint.bin", "r");
	if(Fp != NULL)
	{
  view_count = 0;
  for (int i = 0; i < m_Sinogram->N_theta; i++)
  {
	 
      fread(buffer, sizeof(double), 1, Fp);
	 
  //  if(m_Inputs->ViewMask[i] == 1)
//	  m_Sinogram->InitialGain[view_count++] = (DATA_TYPE)(*buffer);
	  m_Sinogram->InitialGain[i] = (DATA_TYPE)(*buffer);
	  std::cout<<m_Sinogram->InitialGain[i]<<std::endl;
  }
  view_count = 0;
  for (int i = 0; i < m_Sinogram->N_theta; i++)
  {
    fread(buffer, sizeof(double), 1, Fp);
    //if(m_Inputs->ViewMask[i] == 1)
	//m_Sinogram->InitialOffset[view_count++] = (DATA_TYPE)(*buffer);
     m_Sinogram->InitialOffset[i] = (DATA_TYPE)(*buffer);
	 std::cout<<m_Sinogram->InitialOffset[i]<<std::endl;
  }
	}
  fclose(Fp);

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the Gains and Offsets Input file", 0, UpdateProgressMessage);

}
