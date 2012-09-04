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
#include "GainsOffsetsReader.h"
#include "ReconstructionCoreLib/Common/allocate.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsReader::GainsOffsetsReader()
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
  notify("GainsOffsetsReader Starting", 0, UpdateProgressMessage);
  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr inputs = getTomoInputs();


  std::vector<double> fileGains(inputs->fileZSize, inputs->targetGain);//Default value is TARGET_GAIN
  std::vector<double> fileOffsets(inputs->fileZSize, 0);//Defaulted to zero
  std::vector<double> fileVariance(inputs->fileZSize, 1);//Default set to 1

  FILE* Fp = NULL;
  size_t elementsRead = 0;


	//Gains Read
   Fp = fopen(inputs->gainsInputFile.c_str(), "rb");
	if(Fp != NULL)
	{
		// Read all the gains in a single shot
		elementsRead = fread( &(fileGains.front()), sizeof(double), inputs->goodViews.size(), Fp);
		if (elementsRead != inputs->goodViews.size())
		{
			std::cout<<elementsRead<<std::endl;
			std::cout<<inputs->fileZSize<<std::endl;
			setErrorCondition(-1);
			setErrorMessage("Error Reading Gains from File");
			notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
			fclose(Fp);
			return;
		}

		fclose(Fp);
		// Allocate the proper amount of memory for the gains and offsets
		//The normalization and offset parameters for the views
		size_t dims[1] = {sinogram->N_theta};
		sinogram->InitialGain = RealArrayType::New(dims, "InitialGain");

		// Copy just the values of the gains we need from the data read
		// from the file. The indices into the fileGains array are stored
		// in the inputs->goodViews vector

		/*for (size_t i = 0; i < inputs->goodViews.size(); i++)
		{
			sinogram->InitialGain->d[i] = fileGains[inputs->goodViews[i]];
		}*/

		//This has been changed since the # of gains = # of good views
		for (size_t i = 0; i < inputs->goodViews.size(); i++)
		{
			sinogram->InitialGain->d[i] = fileGains[i];
		}

#if 1
		std::cout << "------------Initial Gains-----------" << std::endl;
		for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
		{
			std::cout << "Tilt: " << i_theta<< "  Gain: "<< sinogram->InitialGain->d[i_theta]<< std::endl;
		}

#endif
		setErrorCondition(0);
		setErrorMessage("");
		notify("Done Reading the Gains Input file", 0, UpdateProgressMessage);
	}
	else
	{
		//If no gains file is given just use default values
		std::cout<<"**************************"<<std::endl;
		std::cout<<"Setting Gains to the default value="<<inputs->targetGain<<std::endl;
		std::cout<<"**************************"<<std::endl;
		size_t dims[1] = {sinogram->N_theta};
		sinogram->InitialGain = RealArrayType::New(dims, "InitialGain");
		for(uint16_t i_theta=0; i_theta < inputs->goodViews.size(); i_theta++)
			sinogram->InitialGain->d[i_theta]=inputs->targetGain;
		/*setErrorCondition(-1);
		setErrorMessage("Could not open Gains File");
		notify(getErrorMessage().c_str(), 100, UpdateErrorMessage);
		return;*/
	}

	//Offset Read
	Fp = fopen(inputs->offsetsInputFile.c_str(), "rb");
	if(Fp != NULL)
	{

		elementsRead = fread( &(fileOffsets.front()), sizeof(double), inputs->goodViews.size(), Fp);
		if (elementsRead != inputs->goodViews.size())
		{
			setErrorCondition(-1);
			setErrorMessage("Error Reading Offsets from File");
			notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
			fclose(Fp);
			return;
		}

		fclose(Fp);
		// Allocate the proper amount of memory for the gains and offsets
		//The normalization and offset parameters for the views
		size_t dims[1] = {sinogram->N_theta};
		sinogram->InitialOffset = RealArrayType::New(dims, "InitialOffset");

		// Copy just the values of the gains and offsets we need from the data read
		// from the file. The indices into the fileGains/fileOffsets array are stored
		// in the inputs->goodViews vector
	/*	for (size_t i = 0; i < inputs->goodViews.size(); i++)
		{
			sinogram->InitialOffset->d[i] = fileOffsets[inputs->goodViews[i]];
		}*/

		//This has been changed since the # of gains = # of good views
		for (size_t i = 0; i < inputs->goodViews.size(); i++)
	    {
	    sinogram->InitialOffset->d[i] = fileOffsets[i];
	    }

#if 1

		std::cout << "------------Initial Offsets-----------" << std::endl;
		for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
		{
			std::cout << "Tilt: " << i_theta << "  Offset: "<< sinogram->InitialOffset->d[i_theta] << std::endl;
		}
#endif
		setErrorCondition(0);
		setErrorMessage("");
		notify("Done Reading the Offsets Input file", 0, UpdateProgressMessage);
	}
	else
	{
		setErrorCondition(-1);
		setErrorMessage("Could not open Offsets File");
		notify(getErrorMessage().c_str(), 100, UpdateErrorMessage);
		return;
	}


	//Variance Read
	Fp = fopen(inputs->varianceInputFile.c_str(), "rb");
	if(Fp != NULL)
	{
		// Read all the gains in a single shot
		elementsRead = fread( &(fileVariance.front()), sizeof(double), inputs->goodViews.size(), Fp);
		if (elementsRead != inputs->goodViews.size())
		{
			setErrorCondition(-1);
			setErrorMessage("Error Reading Variance from File");
			notify(getErrorMessage().c_str(), 0, UpdateErrorMessage);
			fclose(Fp);
			return;
		}


		fclose(Fp);
		// Allocate the proper amount of memory for the gains and offsets
		//The normalization and offset parameters for the views
		size_t dims[1] = {sinogram->N_theta};
		sinogram->InitialVariance = RealArrayType::New(dims, "InitialVariance");

		// Copy just the values of the gains and offsets we need from the data read
		// from the file. The indices into the fileGains/fileOffsets array are stored
		// in the inputs->goodViews vector
		/*for (size_t i = 0; i < inputs->goodViews.size(); i++)
		{
			sinogram->InitialVariance->d[i] = fileVariance[inputs->goodViews[i]];
		}*/

		for (size_t i = 0; i < inputs->goodViews.size(); i++)
		{
			sinogram->InitialVariance->d[i] = fileVariance[i];
		}

#if 1
		std::cout << "------------Initial Variance-----------" << std::endl;
		for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
		{
			std::cout << "Tilt: " << i_theta<< "  Variance: "<< sinogram->InitialVariance->d[i_theta]<< std::endl;
		}
#endif
		setErrorCondition(0);
		setErrorMessage("");
		notify("Done Reading the Variance Input file", 0, UpdateProgressMessage);
	}
	else
	{
		std::cout<<"**************************"<<std::endl;
		std::cout<<"Setting Variance to the default value=1"<<std::endl;
		std::cout<<"**************************"<<std::endl;
		size_t dims[1] = {sinogram->N_theta};
		sinogram->InitialVariance = RealArrayType::New(dims, "InitialVariance");
		for(uint16_t i_theta=0; i_theta < inputs->goodViews.size(); i_theta++)
			sinogram->InitialVariance->d[i_theta]=1;
		/*
		setErrorCondition(-1);
		setErrorMessage("Could not open Variance File");
		notify(getErrorMessage().c_str(), 100, UpdateErrorMessage);
		return;*/
	}

}
