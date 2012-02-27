/*
 * ComputeGainsOffets.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: mjackson
 */

#include "ComputeGainsOffsets.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeGainsOffsets::ComputeGainsOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeGainsOffsets::~ComputeGainsOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeGainsOffsets::execute()
{
  // If an error occurs, clean up any memory, call "setErrorCondition(-1)" and
  // also setErrorMessage("Something went wrong"); and then return


	notify("GainsOffsetsCalculation Starting", 0, UpdateProgressMessage);
	SinogramPtr sinogram = getSinogram();//This I assume some how gets the sinogram as it stands now
//	TomoInputs* inputs = getInputs();//This gets the input files

	//The normalization and offset parameters for the views

	size_t dims[3] = {sinogram->N_theta, 0, 0};
	sinogram->InitialGain = RealArrayType::New(dims);
	sinogram->InitialGain->setName("sinogram->InitialGain");
    sinogram->InitialOffset = RealArrayType::New(dims);
    sinogram->InitialOffset->setName("sinogram->InitialOffset");
	sinogram->InitialVariance = RealArrayType::New(dims);
	sinogram->InitialVariance->setName("sinogram->InitialVariance");
	//Form the average Gain per view
	std::vector<DATA_TYPE> AverageGain(sinogram->N_theta,0);
	std::vector<DATA_TYPE> TargetGain(sinogram->N_theta,0);
//	std::vector<std::vector<DATA_TYPE>>LS_Matrix(sinogram->N_theta, std::vector<DATA_TYPE> (2));
	DATA_TYPE** LS_Matrix=(DATA_TYPE**)get_img(2, sinogram->N_theta, sizeof(DATA_TYPE));
	std::vector<DATA_TYPE> LS_Estimates(2,0);

	for(uint16_t i_theta=0; i_theta < sinogram->N_theta; i_theta++)
	{
		DATA_TYPE sum=0;
		for (uint16_t i_r=0; i_r < sinogram->N_r; i_r++) {
			for(uint16_t i_t = 0; i_t <sinogram->N_t; i_t++)
				sum += sinogram->counts->d[i_theta][i_r][i_t];
		}
		sum/=(sinogram->N_r*sinogram->N_t);
		AverageGain[i_theta] = sum;
		TargetGain[i_theta] = 1.0/cos(sinogram->angles[i_theta]*M_PI/180);//Set to 1/cos(tilt_angle)
		LS_Matrix[i_theta][0]=TargetGain[i_theta];
		LS_Matrix[i_theta][1]=1;
	}

	std::cout<<"Average Gains"<<std::endl;
	for(uint16_t i_theta=0; i_theta < sinogram->N_theta; i_theta++)
	{
		std::cout<<AverageGain[i_theta]<<std::endl;
	}

	std::cout<<"Target Gains"<<std::endl;
	for(uint16_t i_theta=0; i_theta < sinogram->N_theta; i_theta++)
	{
		std::cout<<TargetGain[i_theta]<<std::endl;
	}

//	DATA_TYPE max_elt=*max_element(AverageGain.begin(), AverageGain.end());
//	DATA_TYPE min_elt=*min_element(AverageGain.begin(), AverageGain.end());

		//Compute A^t * A
	DATA_TYPE ProdMatrix[2][2];

	for(uint8_t i=0; i < 2;i++)
	for (uint8_t j=0; j < 2; j++)
	{
		ProdMatrix[i][j]=0;
		for (uint8_t k=0; k < sinogram->N_theta; k++) {
			ProdMatrix[i][j]+=(LS_Matrix[k][i]*LS_Matrix[k][j]);
		}
	}



	//Compute T1=(A^t*A)^-1
	DATA_TYPE determinant = ProdMatrix[0][0]*ProdMatrix[1][1] - (ProdMatrix[0][1]*ProdMatrix[0][1]);
	DATA_TYPE InverseProdMatrix[2][2];

	InverseProdMatrix[0][0]=ProdMatrix[1][1]/determinant;
	InverseProdMatrix[1][1]=ProdMatrix[0][0]/determinant;
	InverseProdMatrix[0][1]=-ProdMatrix[1][0]/determinant;
	InverseProdMatrix[1][0]=-ProdMatrix[0][1]/determinant;

	//Compute T2=A^t*Y
	std::vector<DATA_TYPE> TempVector(2,0);
	for (uint8_t j = 0; j < 2; j++) {
	for(uint8_t i_theta=0; i_theta < sinogram->N_theta; i_theta++)
	{
		TempVector[j]+=LS_Matrix[i_theta][j]*AverageGain[i_theta];
		}
	}


	//Compute T1*T2 (2 X 2 matrix with 2 X 1 vector)

	for(uint8_t i=0; i < 2;i++)
	{
	 for (uint8_t k=0; k < 2; k++)
	 {
		 LS_Estimates[i] += InverseProdMatrix[i][k]*TempVector[k];
	 }

	}

	sinogram->TargetGain=TARGET_GAIN;
	std::cout<<"Target Gain"<<sinogram->TargetGain<<std::endl;
	//In this the Gains are all set to the target Gain
	std::cout << "------------Initial Gains-----------" << std::endl;
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	{
		sinogram->InitialGain->d[i_theta] = sinogram->TargetGain;
		std::cout << "Tilt: " << i_theta<< "  Gain: "<< sinogram->InitialGain->d[i_theta]<< std::endl;
	}
	std::cout << "------------Initial Offsets-----------" << std::endl;
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	{
		sinogram->InitialOffset->d[i_theta] = LS_Estimates[1];
		std::cout << "Tilt: " << i_theta << "  Offset: "<< sinogram->InitialOffset->d[i_theta] << std::endl;
	}
	std::cout << "------------Initial Variance-----------" << std::endl;
	for (uint16_t i_theta = 0; i_theta < sinogram->N_theta; i_theta++)
	{
		sinogram->InitialVariance->d[i_theta] = 1;
		std::cout << "Tilt: " << i_theta << "  Variance: "<< sinogram->InitialVariance->d[i_theta] << std::endl;
	}

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done ComputeGainsOffsets", 0, UpdateProgressMessage);
}
