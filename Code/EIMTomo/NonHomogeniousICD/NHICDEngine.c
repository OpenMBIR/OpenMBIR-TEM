/*
 *  ComputationEngineOptimized.c
 *  HAADFSTEM_Xcode
 *
 *  Created by Singanallur Venkatakrishnan on 6/29/11.
 *  Copyright 2011 Purdue University. All rights reserved.
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "EIMTomo/common/EIMTomoTypes.h"
#include "EIMTomo/common/allocate.h"
#include "EIMTomo/common/EMTime.h"
#include "EIMTomo/common/EMMath.h"
#include "NHICDEngine.h"

static char CE_Cancel = 0;



// void CE_cancel()
// {
// CE_Cancel = 1;
// }
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

int CE_MAPICDReconstruct(Sino* Sinogram, Geom* Geometry,CommandLineInputs* CmdInputs)
{

	uint8_t err = 0;
	int16_t Iter,NonHomIter;
	int16_t i,j,k,r,row,col,slice,RowIndex,ColIndex,SliceIndex,Idx,i_r,i_theta,i_t,s,w;
	int16_t NumOfXPixels;
	int32_t q,p;
	//Random Indexing Parameters
	int32_t Index,ArraySize,j_new,k_new;
	int32_t* Counter;
	uint32_t UpdateCount=0;



  //Allocate space for storing columns the A-matrix; an array of pointers to columns
  //AMatrixCol** AMatrix=(AMatrixCol **)get_spc(Geometry->N_x*Geometry->N_z,sizeof(AMatrixCol*));
  //	DetectorResponse = CE_DetectorResponse(0,0,Sinogram,Geometry,VoxelProfile);//System response

#ifdef STORE_A_MATRIX

  AMatrixCol**** AMatrix = (AMatrixCol ****)multialloc(sizeof(AMatrixCol*),3,Geometry->N_y,Geometry->N_z,Geometry->N_x);
#else
  double y;
  double t,tmin,tmax,ProfileThickness;
  int16_t slice_index_min;
  int16_t slice_index_max;
  AMatrixCol*** TempCol = (AMatrixCol***)multialloc(sizeof(AMatrixCol*),2,Geometry->N_z,Geometry->N_x);
  AMatrixCol* Aj;
  AMatrixCol* TempMemBlock;
#endif


	double checksum = 0,temp;
	double RandomIndexj,RandomIndexi;
	double *cost;
	double** VoxelProfile,***DetectorResponse;
	double ***Y_Est;//Estimted Sinogram
	double ***ErrorSino;//Error Sinogram
	double ***Weight;//This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation
    double **UpdateMap;
	double MaxUpdateMagnitude;
	double HammingWindow[5][5]={{0.0013 ,   0.0086 ,   0.0159 ,   0.0086,    0.0013},
		{0.0086 ,   0.0581    ,0.1076   , 0.0581  ,  0.0086},
		{0.0159  ,  0.1076 ,   0.1993 ,   0.1076  ,  0.0159},
		{0.0086  ,  0.0581  ,  0.1076  ,  0.0581  ,  0.0086},
		{0.0013   , 0.0086 ,   0.0159 ,   0.0086  ,  0.0013}};

	FILE *Fp = fopen("ReconstructedSino.bin","w");//Reconstructed Sinogram from initial est

	FILE* Fp2;//Cost function
	FILE *Fp3;//File to store intermediate outputs of reconstruction
	FILE *Fp4 = fopen("UpdateMap.bin","w");

	//Optimization variables
	double low,high;
	double UpdatedVoxelValue;
	double accuracy =1e-7;
	int errorcode;
	//Pointer to  1-D minimization Function
#ifdef WRITE_INTERMEDIATE_RESULTS
	double Fraction = 0.1;//write this fraction of the iterations
	int16_t NumOfWrites = floor((double)(CmdInputs->NumIter)*Fraction);
	int16_t WriteCount = 0;
	char Filename[100];
	char buffer[20];
	void* TempPointer;
	size_t NumOfBytesWritten;
#endif

#ifdef ROI
	uint8_t** Mask;
	double EllipseA,EllipseB,x,z;
#endif

#ifdef COST_CALCULATE
	cost=(double*)get_spc(CmdInputs->NumIter+1,sizeof(double));
#endif

	Y_Est=(double ***)get_3D(Sinogram->N_theta,Sinogram->N_r,Sinogram->N_t,sizeof(double));
	ErrorSino=(double ***)get_3D(Sinogram->N_theta,Sinogram->N_r,Sinogram->N_t,sizeof(double));
	Weight=(double ***)get_3D(Sinogram->N_theta,Sinogram->N_r,Sinogram->N_t,sizeof(double));


	UpdateMap = (double**)get_img(Geometry->N_x, Geometry->N_z, sizeof(double));//width = N_x height = N_z

	for( j = 0; j < Geometry->N_z; j++)
		for (k = 0; k < Geometry->N_x; k++)
		UpdateMap[j][k]=0;


	//calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram Structure are converted to radians
	VoxelProfile = CE_CalculateVoxelProfile(Sinogram,Geometry); //Verified with ML
	CE_CalculateSinCos(Sinogram);
	//Initialize the e-beam
	CE_InitializeBeamProfile(Sinogram); //verified with ML

	//calculate sine and cosine of all angles and store in the global arrays sine and cosine
	DetectorResponse = CE_DetectorResponse(0,0,Sinogram,Geometry,VoxelProfile);//System response

#ifdef ROI
	Mask = (uint8_t **)get_img(Geometry->N_x, Geometry->N_z,sizeof(uint8_t));//width,height
	EllipseA = Geometry->LengthX/2;
	EllipseB = Geometry->LengthZ/2;
	for (i = 0; i < Geometry->N_z ; i++)
	{
		for (j=0; j< Geometry->N_x; j++)
		{
			x = Geometry->x0 + ((double)j + 0.5)*Geometry->delta_xz;
			z = Geometry->z0 + ((double)i+0.5)*Geometry->delta_xz;
			if (((x*x)/(EllipseA*EllipseA)) + ((z*z)/(EllipseB*EllipseB)) < 0.80)
			{
				Mask[i][j] = 1;
			}
			else
			{
				Mask[i][j] =0;
				//for(k = 0;k < Geometry->N_y; k++)
				//Geometry->Object[k][i][j] = 0;
			}

		}
	}


#endif

#ifdef CORRECTION
	//Calculate Normalization factors
	NORMALIZATION_FACTOR = (double*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(double));

	for(i = 0;i < Sinogram->N_r*Sinogram->N_theta;i++)
		NORMALIZATION_FACTOR[i] = 1.0;

#endif

	//Assign the scaling to a global variable
#ifdef BEAM_CALCULATION
	BEAM_WIDTH = (0.5)*Sinogram->delta_r;
#else
	BEAM_WIDTH = Sinogram->delta_r;
#endif

	MRF_P = CmdInputs->p;
	SIGMA_X_P = pow(CmdInputs->SigmaX,MRF_P);

	for(i=0;i<Sinogram->N_theta;i++)
		for(j=0;j< PROFILE_RESOLUTION;j++)
			checksum+=VoxelProfile[i][j];
    printf("CHK SUM%lf\n",checksum);

    checksum=0;


	//Calculating A-Matrix one column at a time
	//For each entry the idea is to initially allocate space for Sinogram.N_theta * Sinogram.N_x
	// And then store only the non zero entries by allocating a new array of the desired size
	//k=0;



	checksum = 0;
	q = 0;

#ifdef STORE_A_MATRIX
	for(i = 0;i < Geometry->N_y; i++)
		for(j = 0;j < Geometry->N_z; j++)
			for (k = 0; k < Geometry->N_x; k++)
			{
				//  AMatrix[q++]=CE_CalculateAMatrixColumn(i,j,Sinogram,Geometry,VoxelProfile);
				AMatrix[i][j][k]=CE_CalculateAMatrixColumn(j,k,i,Sinogram,Geometry,VoxelProfile);//row,col,slice

				for(p = 0; p < AMatrix[i][j][k]->count; p++)
					checksum += AMatrix[i][j][k]->values[p];
				//   printf("(%d,%d,%d) %lf \n",i,j,k,AMatrix[i][j][k]->values);
				checksum = 0;
			}
	printf("Stored A matrix\n");
#else
    for(j=0; j < Geometry->N_z; j++)
    for(k=0; k < Geometry->N_x; k++)
    {
   TempCol[j][k] = CE_CalculateAMatrixColumnPartial(j,k,0,Sinogram,Geometry,DetectorResponse);
//		TempCol[j][k] = CE_CalculateAMatrixColumnPartial(j,k,Sinogram,Geometry,VoxelProfile);
//	printf("%d\n",TempCol[j][k]->count);
    }
#endif

	printf("Geometry-Z %d\n",Geometry->N_z);

	//Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
	//is Y_Est initailized to zero?
	for(i = 0; i < Sinogram->N_theta; i++)
		for(j = 0; j < Sinogram->N_r; j++)
			for(k = 0;k < Sinogram->N_t; k++)
				Y_Est[i][j][k]=0.0;

	/*

	int i,index,ArraySize=400;
	int* Counter = (int*)malloc(400*sizeof(int));
	for(i=0;i < 400;i++)
		Counter[i] = i;
	srand(time(NULL));

	for (i = 0; i < 400; i++)
	{
		index=rand()%ArraySize;
		printf("%d\n",Counter[index]);
		memmove(Counter+index,Counter+index+1,(ArraySize - index)*sizeof(int));
		ArraySize--;
	}*/

		//  p=0;//This is used to access the A-Matrix column correspoding the (i,j,k)the voxel location
	srand(time(NULL));
	ArraySize= Geometry->N_z*Geometry->N_x;
	//ArraySizeK = Geometry->N_x;

	Counter = (int32_t*)malloc(ArraySize*sizeof(int32_t));
	//CounterK = (int*)malloc(ArraySizeK*sizeof(int));

	for(j_new = 0;j_new < ArraySize; j_new++)
		Counter[j_new]=j_new;




//	for(k = 0;k < ArraySizeK; k++)
//		CounterK[k]=k;


		for(j = 0;j < Geometry->N_z; j++)
			for(k = 0;k < Geometry->N_x; k++)
			{
				//generating a random index

				Index = rand()%ArraySize;
				k_new = Counter[Index]%Geometry->N_x;
				j_new = Counter[Index]/Geometry->N_x;
				memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index));
				ArraySize--;

			//	printf("%d,%d\n",j_new,k_new);


			for (i = 0;i < Geometry->N_y; i++)//slice index
				{
#ifdef STORE_A_MATRIX  //A matrix call
				for(q = 0;q < AMatrix[i][j][k]->count; q++)
					//for(q = 0; q < AMatrix[j][k]->count; q++)
				{
					//Convert AMatrix->index[p] to a (slice,row,col)
					row = (int16_t)floor(AMatrix[i][j][k]->index[q]/(Sinogram->N_r*Sinogram->N_t));
					slice =(int16_t) floor((AMatrix[i][j][k]->index[q]%(Sinogram->N_r*Sinogram->N_t))/(Sinogram->N_r));
					col =  (int16_t)(AMatrix[i][j][k]->index[q]%(Sinogram->N_r*Sinogram->N_t))%(Sinogram->N_r);
					//row = (uint32_t)((AMatrix[j][k]->index[q])/Sinogram->N_r);
					//col = (AMatrix[j][k]->index[q])%Sinogram->N_r;
					//printf("(%d,%d)\n",row,col);
					Y_Est[slice][row][col]+=(AMatrix[i][j][k]->values[q] * Geometry->Object[i][j][k]);
					//Y_Est[i][row][col]+=(AMatrix[p]->values[q] * Geometry->Object[i][j][k]);
				}
				//  p++;
#else
                    y = ((double)i+0.5)*Geometry->delta_xy + Geometry->y0;
					t = y;
					tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
				    tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax? t + Geometry->delta_xy/2 : Sinogram->TMax;

		           slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
		           slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

		           if(slice_index_min < 0)
			       slice_index_min = 0;
		           if(slice_index_max >= Sinogram->N_t)
			       slice_index_max = Sinogram->N_t-1;

				for(q = 0;q < TempCol[j_new][k_new]->count; q++)
				{
				    //calculating the footprint of the voxel in the t-direction

					i_theta = floor(TempCol[j_new][k_new]->index[q]/(Sinogram->N_r));
					i_r =  (TempCol[j_new][k_new]->index[q]%(Sinogram->N_r));
					for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
					{
					if(i_t == slice_index_min)
						ProfileThickness = (((i_t+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
					else
					{
						if (i_t == slice_index_max)
						{
							ProfileThickness = (tmax - ((i_t)*Sinogram->delta_t + Sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
						}
						else
						{
							ProfileThickness = Sinogram->delta_t;
						}

					}




					Y_Est[i_theta][i_r][i_t] += (TempCol[j_new][k_new]->values[q] *ProfileThickness* Geometry->Object[j_new][k_new][i]);



					}

				}

#endif
			}

	       }




	//Calculate Error Sinogram - Can this be combined with previous loop?
	//Also compute weights of the diagonal covariance matrix

	for (i_theta=0;i_theta<Sinogram->N_theta;i_theta++)//slice index
	{
		checksum=0;
		for(i_r = 0; i_r < Sinogram->N_r;i_r++)
			for(i_t = 0;i_t < Sinogram->N_t;i_t++)
			{

				// checksum+=Y_Est[i][j][k];
				ErrorSino[i_theta][i_r][i_t] = Sinogram->counts[i_theta][i_r][i_t] - Y_Est[i_theta][i_r][i_t];
				if(Sinogram->counts[i_theta][i_r][i_t] != 0)
					Weight[i_theta][i_r][i_t]=1.0/Sinogram->counts[i_theta][i_r][i_t];
				else
					Weight[i_theta][i_r][i_t]=0;
//#ifdef DEBUG
				//writing the error sinogram

				fwrite(&Y_Est[i_theta][i_r][i_t],sizeof(double),1,Fp);
//#endif
				checksum+=Weight[i_theta][i_r][i_t];
			}
		printf("Check sum of Diagonal Covariance Matrix= %lf\n",checksum);

	}

	free(Y_Est);

	fclose(Fp);





#ifdef COST_CALCULATE
	Fp2 = fopen("CostFunc.bin","w");

	printf("Cost Function Calcluation \n");
	/*********************Initial Cost Calculation***************************************************/
	temp=0;
	cost[0]=0;
	for (i_theta = 0; i_theta < Sinogram->N_theta; i_theta++)
		for (i_r = 0; i_r < Sinogram->N_r; i_r++)
			for( i_t = 0; i_t < Sinogram->N_t; i_t++)
				cost[0] += (ErrorSino[i_theta][i_r][i_t] * ErrorSino[i_theta][i_r][i_t] * Weight[i_theta][i_r][i_t]);
	cost[0]/=2;//accounts for the half in from: (1/2)||Y-AX||^2 + ...

	//Each neighboring pair should be counted only once

	for (i = 0; i < Geometry->N_z; i++)
		for (j = 0; j < Geometry->N_x; j++)
			for(k = 0; k < Geometry->N_y; k++)
			{

				if(k+1 <  Geometry->N_y)
					temp += FILTER[2][1][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j][k+1]),MRF_P);


				if(j+1 < Geometry->N_x)
				{
					if(k-1 >= 0)
						temp += FILTER[0][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k-1]),MRF_P); //<y,z,x>=><z,x,y>


					temp += FILTER[1][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k]),MRF_P);


					if(k+1 < Geometry->N_x)
						temp += FILTER[2][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k+1]),MRF_P);

				}

				if(i+1 < Geometry->N_y)
				{

					if(j-1 >= 0)
						temp += FILTER[1][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k]),MRF_P);

					temp += FILTER[1][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k]),MRF_P);

					if(j+1 < Geometry->N_z)
						temp += FILTER[1][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k]),MRF_P);


					if(j-1 >= 0)
					{
						if(k-1 >= 0)
							temp += FILTER[0][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k-1]),MRF_P);

						if(k+1 < Geometry->N_x)
							temp += FILTER[2][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k+1]),MRF_P);

					}

					if(k-1 >= 0)
						temp += FILTER[0][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k-1]),MRF_P);

					if(j+1 < Geometry->N_z)
					{
						if(k-1 >= 0)
							temp += FILTER[0][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k-1]),MRF_P);

						if(k+1 < Geometry->N_x)
							temp+= FILTER[2][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k+1]),MRF_P);
					}

					if(k+1 < Geometry->N_x)
						temp+= FILTER[2][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k+1]),MRF_P);
				}
			}

	cost[0] += (temp/(MRF_P*SIGMA_X_P));

	printf("%lf\n",cost[0]);

	fwrite(&cost[0],sizeof(double),1,Fp2);
	/*******************************************************************************/
#endif //Cost calculation endif


	//Loop through every voxel updating it by solving a cost function
	srand(time(NULL));

	for(Iter = 0;Iter < CmdInputs->NumIter;Iter++)
	{
		ArraySize = Geometry->N_x*Geometry->N_z;
		for(j_new = 0;j_new < ArraySize; j_new++)
			Counter[j_new]=j_new;

		//printf("Iter %d\n",Iter);

		if(Iter%2  == 0)
		{
		//Homogenous Update
		for(j = 0; j < Geometry->N_z; j++)//Row index
			for(k = 0; k < Geometry->N_x ; k++)//Column index
			{

				Index = rand()%ArraySize;
				k_new = Counter[Index]%Geometry->N_x;
				j_new = Counter[Index]/Geometry->N_x;
				memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index));

				ArraySize--;


				TempMemBlock = TempCol[j_new][k_new];
				UpdateMap[j_new][k_new]=0;

				for (i = 0; i < Geometry->N_y; i++)//slice index
				{
#ifdef ROI
					if (Mask[j_new][k_new] == 1)
					{
#endif
						//Neighborhood of (i,j,k) should be initialized to zeros each time

						for(p = 0; p <= 2; p++)
							for(q = 0; q <= 2; q++)
								for (r = 0; r <= 2;r++)
									NEIGHBORHOOD[p][q][r] = 0.0;

						//For a given (i,j,k) store its 26 point neighborhood
						for(p = -1; p <=1; p++)
							for(q = -1; q <= 1; q++)
								for(r = -1; r <= 1;r++)
								{
									if(i+p >= 0 && i+p < Geometry->N_y)
										if(j_new+q >= 0 && j_new+q < Geometry->N_z)
											if(k_new+r >= 0 && k_new+r < Geometry->N_x)
											{
												NEIGHBORHOOD[p+1][q+1][r+1] = Geometry->Object[q+j_new][r+k_new][p+i];
											}
								}
						NEIGHBORHOOD[1][1][1] = 0.0;

#ifdef DEBUG
						if(i==0 && j == 31 && k == 31)
						{
							printf("***************************\n");
							printf("Geom %lf\n",Geometry->Object[i][31][31]);
							for(p = 0; p <= 2; p++)
								for(q = 0; q <= 2; q++)
									for (r = 0; r <= 2;r++)
										printf("%lf\n",NEIGHBORHOOD[p][q][r]);
						}
#endif
						//Compute theta1 and theta2

						V = Geometry->Object[j_new][k_new][i];//Store the present value of the voxel
						THETA1 = 0.0;
						THETA2 = 0.0;

#ifdef STORE_A_MATRIX

						for (q = 0;q < AMatrix[i][j_new][k_new]->count; q++)
						{

							RowIndex = floor(AMatrix[i][j_new][k_new]->index[q]/(Sinogram->N_r*Sinogram->N_t));
							ColIndex = (AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))%(Sinogram->N_r);
							SliceIndex = floor((AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))/(Sinogram->N_r));

							THETA2 += (AMatrix[i][j_new][k_new]->values[q])*(AMatrix[i][j_new][k_new]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
							THETA1 += ErrorSino[SliceIndex][RowIndex][ColIndex]*(AMatrix[i][j_new][k_new]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
						}
#else

						y = ((double)i+0.5)*Geometry->delta_xy + Geometry->y0;
						t = y;
						tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
						tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax? t + Geometry->delta_xy/2 : Sinogram->TMax;

						slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
						slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

						if(slice_index_min < 0)
							slice_index_min = 0;
						if(slice_index_max >= Sinogram->N_t)
							slice_index_max = Sinogram->N_t-1;

						//TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
						for (q = 0;q < TempMemBlock->count; q++)
						{

							i_theta = floor(TempMemBlock->index[q]/(Sinogram->N_r));
							i_r =  (TempMemBlock->index[q]%(Sinogram->N_r));
							 for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
							 {
								 if(i_t == slice_index_min)
									 ProfileThickness = (((i_t+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
								 else
								 {
									 if (i_t == slice_index_max)
									 {
										 ProfileThickness = (tmax - ((i_t)*Sinogram->delta_t + Sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
									 }
									 else
									 {
										 ProfileThickness = Sinogram->delta_t;
									 }

								 }


								 THETA2 += (ProfileThickness*ProfileThickness)*(TempMemBlock->values[q])*(TempMemBlock->values[q])*Weight[i_theta][i_r][i_t];
								 THETA1 += ErrorSino[i_theta][i_r][i_t]*(TempMemBlock->values[q])*(ProfileThickness)*Weight[i_theta][i_r][i_t];
							 }
						}
#endif

						THETA1*=-1;
						CE_MinMax(&low,&high);


#ifdef DEBUG
						if(i ==0 && j==31 && k==31)
							printf("(%lf,%lf,%lf) \n",low,high,V - (THETA1/THETA2));
#endif

						//Solve the 1-D optimization problem
						//printf("V before updating %lf",V);
						UpdatedVoxelValue = solve(CE_DerivOfCostFunc,low,high,accuracy,&errorcode);
						if(errorcode == 0)
						{
							//    printf("(%lf,%lf,%lf)\n",low,high,UpdatedVoxelValue);
							//printf("%lf\n",UpdatedVoxelValue);
							if(UpdatedVoxelValue < 0.0)//Enforcing positivity constraints
								UpdatedVoxelValue = 0.0;
						}
						else {
							printf("Error \n");
			 	 		}

						//TODO Print appropriate error messages for other values of error code
						Geometry->Object[j_new][k_new][i] = UpdatedVoxelValue;
						//Update the ErrorSinogram

						//Updating the magnitude map for NH-ICD
						UpdateMap[j_new][k_new] +=  fabs(Geometry->Object[j_new][k_new][i] - V);

#ifdef STORE_A_MATRIX
						for (q = 0;q < AMatrix[i][j_new][k_new]->count; q++)
						{

							RowIndex = floor(AMatrix[i][j_new][k_new]->index[q]/(Sinogram->N_r*Sinogram->N_t));
							ColIndex = (AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))%(Sinogram->N_r);
							SliceIndex = floor((AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))/(Sinogram->N_r));
							ErrorSino[SliceIndex][RowIndex][ColIndex] -= (AMatrix[i][j_new][k_new]->values[q] * (Geometry->Object[i][j_new][k_new] - V));
						}
#else

					for (q = 0;q < TempMemBlock->count; q++)
					{

					i_theta = floor(TempMemBlock->index[q]/(Sinogram->N_r));
					i_r =  (TempMemBlock->index[q]%(Sinogram->N_r));
					for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
					{
					if(i_t == slice_index_min)
						ProfileThickness = (((i_t+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//(Sinogram->delta_t); //Will be < Sinogram->delta_t
					else
					{
						if (i_t == slice_index_max)
						{
							ProfileThickness = (tmax - ((i_t)*Sinogram->delta_t + Sinogram->T0));//(Sinogram->delta_t);//Will be < Sinogram->delta_t
						}
						else
						{
							ProfileThickness = (Sinogram->delta_t);
						}

					}

							ErrorSino[i_theta][i_r][i_t] -= (TempMemBlock->values[q] *ProfileThickness* (Geometry->Object[j_new][k_new][i] - V));
					}
					}
#endif
						Idx++;



#ifdef ROI
				} //closing the if for region of interest only updates
#endif

				}
			}
			//Need to filter Update map
			for(j = 0; j < Geometry->N_z; j ++)
				for(k = 0; k < Geometry->N_x; k++)
				{
					temp=0;
					for(s=-2 ; s <= 2 ; s++)
						for(w = -2; w <= 2; w++)
							if(j+s >= 0 && j+s < Geometry->N_z && k + w >= 0 && k+w < Geometry->N_x)
								temp+=(UpdateMap[j+s][k+w]*HammingWindow[s+2][w+2]);
					UpdateMap[j][k]=temp;

				}

			//Find Max
			MaxUpdateMagnitude = -1;
			for (j=20; j < Geometry->N_z-20; j++)
				for (k=0+50; k < Geometry->N_x-50; k++)
					if(MaxUpdateMagnitude < UpdateMap[j][k])
						MaxUpdateMagnitude = UpdateMap[j][k];






		}
		else
		{

			//Non-homogenous Update - K times of the values with highest say 10% of  UpdateMap vals

			for(NonHomIter = 0; NonHomIter < 10;NonHomIter++)//Need to change 10 to some other number possibly
			{
				ArraySize = Geometry->N_x*Geometry->N_z;
				for(j_new = 0;j_new < ArraySize; j_new++)
					Counter[j_new]=j_new;
				//Need to filter the UpdateMap
				UpdateCount=0;

				//printf("MaxUpdateMag=%lf\n",MaxUpdateMagnitude);

			for(j = 0; j < Geometry->N_z; j++)//Row index
				for(k = 0; k < Geometry->N_x ; k++)//Column index
				{
					Index = rand()%ArraySize;
					k_new = Counter[Index]%Geometry->N_x;
					j_new = Counter[Index]/Geometry->N_x;
					memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index));

					ArraySize--;


					TempMemBlock = TempCol[j_new][k_new];

					if(UpdateMap[j_new][k_new] > 0.4*MaxUpdateMagnitude)
					{
					UpdateCount++;
					UpdateMap[j_new][k_new]=0;
					for (i = 0; i < Geometry->N_y; i++)//slice index
					{
#ifdef ROI
						if (Mask[j_new][k_new] == 1)
						{
#endif

							//Neighborhood of (i,j,k) should be initialized to zeros each time

							for(p = 0; p <= 2; p++)
								for(q = 0; q <= 2; q++)
									for (r = 0; r <= 2;r++)
										NEIGHBORHOOD[p][q][r] = 0.0;

							//For a given (i,j,k) store its 26 point neighborhood
							for(p = -1; p <=1; p++)
								for(q = -1; q <= 1; q++)
									for(r = -1; r <= 1;r++)
									{
										if(i+p >= 0 && i+p < Geometry->N_y)
											if(j_new+q >= 0 && j_new+q < Geometry->N_z)
												if(k_new+r >= 0 && k_new+r < Geometry->N_x)
												{
													NEIGHBORHOOD[p+1][q+1][r+1] = Geometry->Object[q+j_new][r+k_new][p+i];
												}
									}
							NEIGHBORHOOD[1][1][1] = 0.0;

#ifdef DEBUG
							if(i==0 && j == 31 && k == 31)
							{
								printf("***************************\n");
								printf("Geom %lf\n",Geometry->Object[i][31][31]);
								for(p = 0; p <= 2; p++)
									for(q = 0; q <= 2; q++)
										for (r = 0; r <= 2;r++)
											printf("%lf\n",NEIGHBORHOOD[p][q][r]);
							}
#endif
							//Compute theta1 and theta2

							V = Geometry->Object[j_new][k_new][i];//Store the present value of the voxel
							THETA1 = 0.0;
							THETA2 = 0.0;

#ifdef STORE_A_MATRIX

							for (q = 0;q < AMatrix[i][j_new][k_new]->count; q++)
							{

								RowIndex = floor(AMatrix[i][j_new][k_new]->index[q]/(Sinogram->N_r*Sinogram->N_t));
								ColIndex = (AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))%(Sinogram->N_r);
								SliceIndex = floor((AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))/(Sinogram->N_r));

								THETA2 += (AMatrix[i][j_new][k_new]->values[q])*(AMatrix[i][j_new][k_new]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
								THETA1 += ErrorSino[SliceIndex][RowIndex][ColIndex]*(AMatrix[i][j_new][k_new]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
							}
#else

							y = ((double)i+0.5)*Geometry->delta_xy + Geometry->y0;
							t = y;
							tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
							tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax? t + Geometry->delta_xy/2 : Sinogram->TMax;

							slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
							slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

							if(slice_index_min < 0)
								slice_index_min = 0;
							if(slice_index_max >= Sinogram->N_t)
								slice_index_max = Sinogram->N_t-1;

							//TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
							for (q = 0;q < TempMemBlock->count; q++)
							{

								i_theta = floor(TempMemBlock->index[q]/(Sinogram->N_r));
								i_r =  (TempMemBlock->index[q]%(Sinogram->N_r));
								for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
								{
									if(i_t == slice_index_min)
										ProfileThickness = (((i_t+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
									else
									{
										if (i_t == slice_index_max)
										{
											ProfileThickness = (tmax - ((i_t)*Sinogram->delta_t + Sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
										}
										else
										{
											ProfileThickness = Sinogram->delta_t;
										}

									}


									THETA2 += (ProfileThickness*ProfileThickness)*(TempMemBlock->values[q])*(TempMemBlock->values[q])*Weight[i_theta][i_r][i_t];
									THETA1 += ErrorSino[i_theta][i_r][i_t]*(TempMemBlock->values[q])*(ProfileThickness)*Weight[i_theta][i_r][i_t];
								}
							}
#endif

							THETA1*=-1;
							CE_MinMax(&low,&high);


#ifdef DEBUG
							if(i ==0 && j==31 && k==31)
								printf("(%lf,%lf,%lf) \n",low,high,V - (THETA1/THETA2));
#endif

							//Solve the 1-D optimization problem
							//printf("V before updating %lf",V);
							UpdatedVoxelValue = solve(CE_DerivOfCostFunc,low,high,accuracy,&errorcode);
							if(errorcode == 0)
							{
								//    printf("(%lf,%lf,%lf)\n",low,high,UpdatedVoxelValue);
								//printf("%lf\n",UpdatedVoxelValue);
								if(UpdatedVoxelValue < 0.0)//Enforcing positivity constraints
									UpdatedVoxelValue = 0.0;
							}
							else {
								printf("Error \n");
							}

							//TODO Print appropriate error messages for other values of error code
							Geometry->Object[j_new][k_new][i] = UpdatedVoxelValue;
							//Update the ErrorSinogram

							//Updating the magnitude map for NH-ICD
							UpdateMap[j_new][k_new] +=  fabs(Geometry->Object[j_new][k_new][i] - V);

#ifdef STORE_A_MATRIX
							for (q = 0;q < AMatrix[i][j_new][k_new]->count; q++)
							{

								RowIndex = floor(AMatrix[i][j_new][k_new]->index[q]/(Sinogram->N_r*Sinogram->N_t));
								ColIndex = (AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))%(Sinogram->N_r);
								SliceIndex = floor((AMatrix[i][j_new][k_new]->index[q]%(Sinogram->N_r*Sinogram->N_t))/(Sinogram->N_r));
								ErrorSino[SliceIndex][RowIndex][ColIndex] -= (AMatrix[i][j_new][k_new]->values[q] * (Geometry->Object[i][j_new][k_new] - V));
							}
#else

							for (q = 0;q < TempMemBlock->count; q++)
							{

								i_theta = floor(TempMemBlock->index[q]/(Sinogram->N_r));
								i_r =  (TempMemBlock->index[q]%(Sinogram->N_r));
								for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
								{
									if(i_t == slice_index_min)
										ProfileThickness = (((i_t+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//(Sinogram->delta_t); //Will be < Sinogram->delta_t
									else
									{
										if (i_t == slice_index_max)
										{
											ProfileThickness = (tmax - ((i_t)*Sinogram->delta_t + Sinogram->T0));//(Sinogram->delta_t);//Will be < Sinogram->delta_t
										}
										else
										{
											ProfileThickness = (Sinogram->delta_t);
										}

									}

									ErrorSino[i_theta][i_r][i_t] -= (TempMemBlock->values[q] *ProfileThickness* (Geometry->Object[j_new][k_new][i] - V));
								}
							}
#endif
							Idx++;


#ifdef ROI
					} //closing the if for region of interest only updates
#endif


				    }
						}
				}
				//printf("Number of voxel lines updated =%d\n",UpdateCount);
				//Need to filter Update map
				for(j = 0; j < Geometry->N_z; j ++)
					for(k = 0; k < Geometry->N_x; k++)
					{
						temp=0;
						for(s=-2 ; s <= 2 ; s++)
							for(w = -2; w <= 2; w++)
								if(j+s >= 0 && j+s < Geometry->N_z && k + w >= 0 && k+w < Geometry->N_x)
									temp+=(UpdateMap[j+s][k+w]*HammingWindow[s+2][w+2]);
						UpdateMap[j][k]=temp;

					}

				MaxUpdateMagnitude = -1;
				for (j=0+20; j < Geometry->N_z-20; j++)
					for (k=0+50; k < Geometry->N_x-50; k++)
						if(UpdateMap[j][k] > MaxUpdateMagnitude)
							MaxUpdateMagnitude = UpdateMap[j][k];





			}

		}

		if (CE_Cancel == 1)
		{
			return err;
		}





#ifdef COST_CALCULATE
		/*********************Cost Calculation***************************************************/
		temp=0;
		cost[Iter+1]=0;
		for (i_theta = 0; i_theta < Sinogram->N_theta; i_theta++)
			for (i_r = 0; i_r < Sinogram->N_r; i_r++)
				for( i_t = 0; i_t < Sinogram->N_t; i_t++)
					cost[Iter+1] += (ErrorSino[i_theta][i_r][i_t] * ErrorSino[i_theta][i_r][i_t] * Weight[i_theta][i_r][i_t]);
		cost[Iter+1]/=2;//Accounting for (1/2) in the cost function
		//Each neighboring pair should be counted only once

		for (i = 0; i < Geometry->N_z; i++)
			for (j = 0; j < Geometry->N_x; j++)
				for(k = 0; k < Geometry->N_y; k++)
				{

					if(k+1 <  Geometry->N_y)
						temp += FILTER[2][1][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j][k+1]),MRF_P);


					if(j+1 < Geometry->N_x)
					{
						if(k-1 >= 0)
							temp += FILTER[0][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k-1]),MRF_P);


						temp += FILTER[1][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k]),MRF_P);


						if(k+1 < Geometry->N_y)
							temp += FILTER[2][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k+1]),MRF_P);

					}

					if(i+1 < Geometry->N_z)
					{

						if(j-1 >= 0)
							temp += FILTER[1][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k]),MRF_P);

						temp += FILTER[1][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k]),MRF_P);

						if(j+1 < Geometry->N_x)
							temp += FILTER[1][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k]),MRF_P);


						if(j-1 >= 0)
						{
							if(k-1 >= 0)
								temp += FILTER[0][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k-1]),MRF_P);

							if(k+1 < Geometry->N_y)
								temp += FILTER[2][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k+1]),MRF_P);

						}

						if(k-1 >= 0)
							temp += FILTER[0][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k-1]),MRF_P);

						if(j+1 < Geometry->N_x)
						{
							if(k-1 >= 0)
								temp += FILTER[0][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k-1]),MRF_P);

							if(k+1 < Geometry->N_y)
								temp+= FILTER[2][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k+1]),MRF_P);
						}

						if(k+1 < Geometry->N_y)
							temp+= FILTER[2][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k+1]),MRF_P);
					}
				}

		cost[Iter] += (temp/(MRF_P*SIGMA_X_P));


		//Calculating the approximate error magnitude in terms of the updates
/*		temp=0;
		for (j=0; j < Geometry->N_z; j++)
			for (k=0; k < Geometry->N_x; k++)
				temp+=UpdateMap[j][k];

		cost[Iter] = temp;*/

		printf("%lf\n",cost[Iter]);

		fwrite(&cost[Iter],sizeof(double),1,Fp2);
		/*******************************************************************************/
#else
		printf("%d\n",Iter);
#endif //Cost calculation endif

#ifdef WRITE_INTERMEDIATE_RESULTS

		if(Iter == NumOfWrites*WriteCount)
		{
			WriteCount++;
			sprintf(buffer,"%d",Iter);
			sprintf(Filename,"ReconstructedObjectAfterIter");
			strcat(Filename,buffer);
			strcat(Filename,".bin");
			Fp3 = fopen(Filename, "w");
			//	for (i=0; i < Geometry->N_y; i++)
			//		for (j=0; j < Geometry->N_z; j++)
			//			for (k=0; k < Geometry->N_x; k++)
			TempPointer = Geometry->Object;
			NumOfBytesWritten=fwrite(&(Geometry->Object[0][0][0]), sizeof(double),Geometry->N_x*Geometry->N_y*Geometry->N_z, Fp3);
	 		printf("%d\n",NumOfBytesWritten);


			fclose(Fp3);
		}
#endif


	}

	for (j=0; j < Geometry->N_z; j++)
		for (k=0; k < Geometry->N_x; k++)
			fwrite(&UpdateMap[j][k], sizeof(double), 1, Fp4);

	free_img((void*)VoxelProfile);
	//free(AMatrix);
#ifdef STORE_A_MATRIX
	multifree(AMatrix,2);
	//#else
	//	free((void*)TempCol);
#endif
	free_3D((void*)ErrorSino);
	free_3D((void*)Weight);
#ifdef COST_CALCULATE
	fclose(Fp2);// writing cost function
#endif
	//free_3D(neighborhood);
	// Get values from ComputationInputs and perform calculation
	// Return any error code
	return err;
}


/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void CE_MinMax(double *low,double *high)
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



	*low = (*low > (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *low);
	*high = (*high < (V - (THETA1/THETA2)) ? (V - (THETA1/THETA2)): *high);
}



void* CE_CalculateVoxelProfile(Sino *Sinogram,Geom *Geometry)
{
	double angle,MaxValLineIntegral;
	double temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;
	double** VoxProfile = (double**)multialloc(sizeof(double),2,Sinogram->N_theta,PROFILE_RESOLUTION);
	double checksum=0;
	uint16_t i,j;
	FILE* Fp = fopen("VoxelProfile.bin","w");

	for (i=0;i<Sinogram->N_theta;i++)
	{
		Sinogram->angles[i]=Sinogram->angles[i]*(PI/180.0);
		angle=Sinogram->angles[i];
		while(angle > PI/2)
			angle -= PI/2;

		while(angle < 0)
			angle +=PI/2;

		if(angle <= PI/4)
		{
			MaxValLineIntegral = Geometry->delta_xz/cos(angle);
		}
		else
		{
			MaxValLineIntegral = Geometry->delta_xz/cos(PI/2-angle);
		}
		temp=cos(PI/4);
		dist1 = temp * cos((PI/4.0 - angle));
		dist2 = temp * fabs((cos((PI/4.0 + angle))));
		LeftCorner = 1-dist1;
		LeftNear = 1-dist2;
		RightNear = 1+dist2;
		RightCorner = 1+dist1;

		for(j = 0;j<PROFILE_RESOLUTION;j++)
		{
			t = 2.0*j / PROFILE_RESOLUTION;//2 is the normalized length of the profile (basically equl to 2*delta_xz)
			if(t <= LeftCorner || t >= RightCorner)
				VoxProfile[i][j] = 0;
			else if(t > RightNear)
				VoxProfile[i][j] = MaxValLineIntegral*(RightCorner-t)/(RightCorner-RightNear);
			else if(t >= LeftNear)
				VoxProfile[i][j] = MaxValLineIntegral;
			else
				VoxProfile[i][j] = MaxValLineIntegral*(t-LeftCorner)/(LeftNear-LeftCorner);

			fwrite(&VoxProfile[i][j],sizeof(double),1,Fp);
			checksum+=VoxProfile[i][j];
		}

	}

	//printf("Pixel Profile Check sum =%lf\n",checksum);
	fclose(Fp);
	return VoxProfile;
}

/*******************************************************************
 Find each column of the Forward Projection Matrix A
 ********************************************************************/
void ForwardProject(Sino *Sinogram,Geom* Geom)
{

}

void* CE_CalculateAMatrixColumn(uint16_t row,uint16_t col, uint16_t slice, Sino* Sinogram,Geom* Geometry,double** VoxelProfile)
{
	int32_t i,j,k,sliceidx;
	double x,z,y;
	double r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
	double t; //this is similar to r but along the y direction
	double tmin,tmax;
	double rmax,rmin;//stores the start and end points of the pixel profile on the detector
	double RTemp,TempConst,checksum = 0,Integral = 0;
	double LeftEndOfBeam;
	double MaximumSpacePerColumn;//we will use this to allocate space
	double AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
	double ProfileThickness;
	int32_t index_min,index_max,slice_index_min,slice_index_max;//stores the detector index in which the profile lies
	int32_t BaseIndex,FinalIndex,ProfileIndex=0;
	uint32_t count = 0;



#ifdef DISTANCE_DRIVEN
	double d1,d2; //These are the values of the detector boundaries
#endif

	AMatrixCol* Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

	// printf("Space allocated for column %d %d\n",row,col);

	//Temp->index = (uint32_t*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(uint32_t));
	//Temp->values = (double*)multialloc(sizeof(double),1,Sinogram->N_r*Sinogram->N_theta);//makes the values =0

	x = Geometry->x0 + ((double)col+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	z = Geometry->z0 + ((double)row+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	y = Geometry->y0 + ((double)slice + 0.5)*Geometry->delta_xy;

	TempConst=(PROFILE_RESOLUTION)/(2*Geometry->delta_xz);


	//	Temp->values = (double*)calloc(Sinogram->N_t*Sinogram->N_r*Sinogram->N_theta,sizeof(double));//(double*)get_spc(Sinogram->N_r*Sinogram->N_theta,sizeof(double));//makes the values =0

	//alternately over estimate the maximum size require for a single AMatrix column
	AvgNumXElements = ceil(3*Geometry->delta_xz/Sinogram->delta_r);
	AvgNumYElements = ceil(3*Geometry->delta_xy/Sinogram->delta_t);
	MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*Sinogram->N_theta;

	Temp->values = (double*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(double));
	Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));

	//printf("%lf",Temp->values[10]);

#ifdef AREA_WEIGHTED
	for(i=0;i<Sinogram->N_theta;i++)
	{

		r = x*cosine[i] - z*sine[i];
		t = y;

		rmin = r - Geometry->delta_xz;
		rmax = r + Geometry->delta_xz;

		tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
		tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax ? t + Geometry->delta_xy/2 : Sinogram->TMax;

		if(rmax < Sinogram->R0 || rmin > Sinogram->RMax)
			continue;



		index_min = floor(((rmin - Sinogram->R0)/Sinogram->delta_r));
		index_max = floor((rmax - Sinogram->R0)/Sinogram->delta_r);

		if(index_max >= Sinogram->N_r)
			index_max = Sinogram->N_r - 1;

		if(index_min < 0)
			index_min = 0;

		slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
		slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= Sinogram->N_t)
			slice_index_max = Sinogram->N_t -1;

		BaseIndex = i*Sinogram->N_r*Sinogram->N_t;

		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);


		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
		//of same dimesions


		for(j = index_min;j <= index_max; j++)//Check
		{

			Integral = 0.0;

			//Accounting for Beam width
			RTemp = (Sinogram->R0 + (((double)j) + 0.5) *(Sinogram->delta_r));//the 0.5 is to get to the center of the detector

			LeftEndOfBeam = RTemp - (BEAM_WIDTH/2);//Consider pre storing the divide by 2

			//   if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
			//printf("%d %lf %lf\n",j, LeftEndOfBeam,rmin);


			for(k=0; k < BEAM_RESOLUTION; k++)
			{

				RTemp = LeftEndOfBeam + ((((double)k)*(BEAM_WIDTH))/BEAM_RESOLUTION);

				if (RTemp-rmin >= 0.0)
				{
					ProfileIndex = (int32_t)floor((RTemp-rmin)*TempConst);//Finding the nearest neighbor profile to the beam
					//if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
					//printf("%d\n",ProfileIndex);
					if(ProfileIndex > PROFILE_RESOLUTION)
						ProfileIndex = PROFILE_RESOLUTION;
				}
				if(ProfileIndex < 0)
					ProfileIndex = 0;


				if(ProfileIndex >= 0 && ProfileIndex < PROFILE_RESOLUTION)
				{
#ifdef BEAM_CALCULATION
					Integral+=(BeamProfile[k]*VoxelProfile[i][ProfileIndex]);
#else
					Integral+=(VoxelProfile[i][ProfileIndex]/PROFILE_RESOLUTION);
#endif
					//	if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
					//	 printf("Index %d %lf Voxel %lf I=%d\n",ProfileIndex,BeamProfile[k],VoxelProfile[2][274],i);
				}

			}
			if(Integral > 0.0)
			{
				//	printf("Entering, Final Index %d %d\n",FinalIndex,Temp->values[0]);
				//		printf("Done %d %d\n",slice_index_min,slice_index_max);
				for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
				{
					if(sliceidx == slice_index_min)
						ProfileThickness = (((sliceidx+1)*Sinogram->delta_t + Sinogram->T0) - tmin);//Sinogram->delta_t; //Will be < Sinogram->delta_t
					else
					{
						if (sliceidx == slice_index_max)
						{
							ProfileThickness = (tmax - ((sliceidx)*Sinogram->delta_t + Sinogram->T0));//Sinogram->delta_t;//Will be < Sinogram->delta_t
						}
						else
						{
							ProfileThickness = Sinogram->delta_t;//Sinogram->delta_t;
						}

					}
					if(ProfileThickness > 0)
					{

						FinalIndex = BaseIndex + (int32_t)j + (int32_t)sliceidx * Sinogram->N_r;
						Temp->values[count] = Integral*ProfileThickness;
						Temp->index[count] = FinalIndex;//can instead store a triple (row,col,slice) for the sinogram
						//printf("Done\n");
#ifdef CORRECTION
						Temp->values[count]/=NORMALIZATION_FACTOR[FinalIndex];//There is a normalizing constant for each measurement . So we have a total of Sinogram.N_x * Sinogram->N_theta values
#endif
						//printf("Access\n");
						count++;
					}

				}
			}

			//  End of Beam width accounting


			//        ProfileIndex=(uint32_t)(((RTemp-rmin)*TempConst));
			// //    //   printf("%d\n",ProfileIndex);
			//        if(ProfileIndex>=0 && ProfileIndex < PROFILE_RESOLUTION)
			//        {
			//  	if(VoxelProfile[i][ProfileIndex] > 0.0)
			//  	{
			//          Temp->values[FinalIndex]=VoxelProfile[i][ProfileIndex];
			//         // Temp->index[count++]=FinalIndex;
			// 	 count++;
			//  	}
			//        }
		}



	}

#endif

#ifdef DISTANCE_DRIVEN

	for(i=0;i<Sinogram->N_theta;i++)
	{

		r = x*cosine[i] - z*sine[i];
		t = y;

		tmin = (t - Geometry->delta_xy/2) > -Geometry->LengthY/2 ? t-Geometry->delta_xy/2 : -Geometry->LengthY/2;
		tmax = (t + Geometry->delta_xy/2) <= Geometry->LengthY/2 ? t + Geometry->delta_xy/2 : Geometry->LengthY/2;


		if(Sinogram->angles[i]*(180/PI) >= -45 && Sinogram->angles[i]*(180/PI) <= 45)
		{
			rmin = r - (Geometry->delta_xz/2)*(cosine[i]);
			rmax = r + (Geometry->delta_xz/2)*(cosine[i]);
		}

		else
		{
			rmin = r - (Geometry->delta_xz/2)*fabs(sine[i]);
			rmax = r + (Geometry->delta_xz/2)*fabs(sine[i]);
		}


		if(rmax < Sinogram->R0 || rmin > Sinogram->R0 + Sinogram->N_r*Sinogram->delta_r)
			continue;

		index_min = floor(((rmin-Sinogram->R0)/Sinogram->delta_r));
		index_max = floor((rmax-Sinogram->R0)/Sinogram->delta_r);

		slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
		slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= Sinogram->N_t)
			slice_index_max = Sinogram->N_t -1;


		if(index_max >= Sinogram->N_r)
			index_max = Sinogram->N_r-1;

		if(index_min < 0)
			index_min = 0;

		BaseIndex = i*Sinogram->N_r*Sinogram->N_t;

		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);



		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
		//of same dimesions

		for(j = index_min;j <= index_max; j++)//Check
		{
			d1  = ((double)j)*Sinogram->delta_r + Sinogram->R0;
			d2 =  d1 + Sinogram->delta_r;

			if(rmax < d1)
			{
				Integral = 0;
			}
			else
			{
				if(rmin > d1 && rmin < d2 && rmax > d2)
				{
					Integral = (d2 - rmin)/Sinogram->delta_r;
				}
				else
				{
					if(rmin >= d1 && rmin <= d2 && rmax >= d1 && rmax <= d2)
						Integral= (rmax - rmin)/Sinogram->delta_r;
					else
						if(rmax > d1 && rmax < d2 && rmin < d1)
							Integral = (rmax - d1)/Sinogram->delta_r;
						else
							if( rmin < d1 && rmax > d2)
								Integral = 1;
							else
								Integral = 0;
				}


			}
			if(Integral > 0)
			{
				//   printf("Final Index %d %lf\n",FinalIndex,cosine[i]);
				for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
				{
					if(sliceidx == slice_index_min)
						ProfileThickness = (((sliceidx+1)*Sinogram->delta_t + Sinogram->T0) - tmin)*Sinogram->delta_t;
					else {
						if (sliceidx == slice_index_max)
						{
							ProfileThickness = (tmax - ((sliceidx)*Sinogram->delta_t + Sinogram->T0))*Sinogram->delta_t;
						}
						else {
							ProfileThickness = Sinogram->delta_t;
						}

					}
					if (ProfileThickness > 0)
					{
						FinalIndex = BaseIndex + (uint32_t)j + (int32_t)sliceidx * Sinogram->N_r;

						Temp->values[count] = Integral*ProfileThickness;
						Temp->index[count] = FinalIndex;
#ifdef CORRECTION
						Temp->values[count]/=NORMALIZATION_FACTOR[FinalIndex];//There is a normalizing constant for each measurement . So we have a total of Sinogram.N_x * Sinogram->N_theta values
#endif
						count++;
					}
				}
			}
			else
			{
				//  printf("%lf \n",Sinogram->angles[i]*180/PI);
			}
		}



	}


#endif
	// printf("Final Space allocation for column %d %d\n",row,col);

	Ai->values=(double*)get_spc(count,sizeof(double));
	Ai->index=(uint32_t*)get_spc(count,sizeof(uint32_t));
	k=0;
	for(i = 0; i < count; i++)
	{
		if(Temp->values[i] > 0.0)
		{
			Ai->values[k]=Temp->values[i];
			checksum+=Ai->values[k];
			Ai->index[k++]=Temp->index[i];
		}

	}

	Ai->count=k;

	//printf("%d %d \n",Ai->count,count);

	//printf("(%d,%d) %lf \n",row,col,checksum);

	free(Temp->values);
	free(Temp->index);
	free(Temp);
	return Ai;

}

/* Initializes the global variables cosine and sine to speed up computation
 */
void CE_CalculateSinCos(Sino* Sinogram)
{
	uint16_t i;
	cosine=(double*)get_spc(Sinogram->N_theta,sizeof(double));
	sine=(double*)get_spc(Sinogram->N_theta,sizeof(double));
	for(i=0;i<Sinogram->N_theta;i++)
	{
		cosine[i]=cos(Sinogram->angles[i]);
		sine[i]=sin(Sinogram->angles[i]);

	}
}

void CE_InitializeBeamProfile(Sino* Sinogram)
{
	uint16_t i;
	double sum=0,W;
	BeamProfile=(double*)get_spc(BEAM_RESOLUTION,sizeof(double));
	W=BEAM_WIDTH/2;
	for (i=0; i < BEAM_RESOLUTION ;i++)
	{
		//BeamProfile[i] = (1.0/(BEAM_WIDTH)) * ( 1 + cos ((PI/W)*fabs(-W + i*(BEAM_WIDTH/BEAM_RESOLUTION))));
		BeamProfile[i] = 0.54 - 0.46*cos((2*PI/BEAM_RESOLUTION)*i);
		sum=sum+BeamProfile[i];
	}

	//Normalize the beam to have an area of 1

	for (i=0; i < BEAM_RESOLUTION ;i++)
	{

		BeamProfile[i]/=sum;
		// printf("%lf\n",BeamProfile[i]);
	}



}

double CE_DerivOfCostFunc(double u)
{
	double temp=0,value;
	uint8_t i,j,k;
	for (i = 0;i < 3;i++)
		for (j = 0; j <3; j++)
			for (k = 0;k < 3; k++)
			{

				if( u - NEIGHBORHOOD[i][j][k] >= 0.0)
					temp+=(FILTER[i][j][k]*(1.0)*pow(fabs(u-NEIGHBORHOOD[i][j][k]),(MRF_P-1)));
				else
					temp+=(FILTER[i][j][k]*(-1.0)*pow(fabs(u-NEIGHBORHOOD[i][j][k]),(MRF_P -1)));
			}

	//printf("V WHile updating %lf\n",V);
	//scanf("Enter value %d\n",&k);
	value = THETA1 +  THETA2*(u-V) + (temp/SIGMA_X_P);

	return value;
}


double  solve(
			  double (*f)(), /* pointer to function to be solved */
			  double a,      /* minimum value of solution */
			  double b,      /* maximum value of solution */
			  double err,    /* accuarcy of solution */
			  int *code      /* error code */
			  )
/* Solves equation (*f)(x) = 0 on x in [a,b]. Uses half interval method.*/
/* Requires that (*f)(a) and (*f)(b) have opposite signs.		*/
/* Returns code=0 if signs are opposite.				*/
/* Returns code=1 if signs are both positive. 				*/
/* Returns code=-1 if signs are both negative. 				*/
{
	int     signa,signb,signc;
	double  fa,fb,fc,c,signaling_nan();
	double  dist;

	fa = (*f)(a);
	signa = fa>0;
	fb = (*f)(b);
	signb = fb>0;

	/* check starting conditions */
	if( signa==signb ) {
		if(signa==1) *code = 1;
		else *code = -1;
		return(0.0);
	}
	else *code = 0;

	/* half interval search */
	if( (dist=b-a)<0 ) dist = -dist;
	while(dist>err) {
		c = (b+a)/2;
		fc = (*f)(c);  signc = fc>0;
		if(signa == signc) { a = c; fa = fc; }
		else { b = c; fb = fc; }
		if( (dist=b-a)<0 ) dist = -dist;
	}

	/* linear interpolation */
	if( (fb-fa)==0 ) return(a);
	else {
		c = (a*fb - b*fa)/(fb-fa);
		return(c);
	}
}

double CE_ComputeCost(double*** ErrorSino,double*** Weight,Sino* Sinogram,Geom* Geometry)
{
	double cost=0,temp=0;
	int16_t i,j,k,p,q,r;

	//printf("Cost calculation..\n");
    //  printf("%lf %lf\n",ErrorSino[5][5][5],Weight[5][5][5]);

	for (i = 0; i < Sinogram->N_t; i++)
		for (j = 0; j < Sinogram->N_theta; j++)
			for( k = 0; k < Sinogram->N_r; k++)
				cost+=(ErrorSino[i][j][k] * ErrorSino[i][j][k] * Weight[i][j][k]);

	for( i=0; i < Geometry->N_y;i++ )
		for (j = 0; j < Geometry->N_z; j++)
			for( k = 0; k < Geometry->N_x; k++)
			{
				for(p = -1; p <= 1; p++)
					for(q = -1; q <= 1; q++)
						for(r = -1; r <= 1; r++)
							if(i+p >= 0 && i+p < Geometry->N_y && j+q >= 0 && j+q < Geometry->N_z && k+r >= 0 && k+r < Geometry->N_x)
								temp+=(FILTER[p+1][q+1][r+1]*pow(fabs(Geometry->Object[i][j][k] - Geometry->Object[i+p][j+q][k+r]),MRF_P));

			}

	//printf("Cost calculation End..\n");

	cost+=(temp/(MRF_P*SIGMA_X_P));
	return cost;
}

void* CE_DetectorResponse(uint16_t row,uint16_t col,Sino* Sinogram,Geom* Geometry,double** VoxelProfile)
{
	FILE* Fp = fopen("DetectorResponse.bin","w");
	double Offset,Temp,r,sum=0,rmin,ProfileCenterR,ProfileCenterT,TempConst,t,tmin;
	double OffsetR,OffsetT,Left;
	double ***H;
	double r0 = -(BEAM_WIDTH)/2;
	double t0 = -(BEAM_WIDTH)/2;
	double StepSize = BEAM_WIDTH/BEAM_RESOLUTION;
	int16_t i,j,k,p,l,NumOfDisplacements,ProfileIndex;
	//NumOfDisplacements=32;
	H = (double***)get_3D(1, Sinogram->N_theta,DETECTOR_RESPONSE_BINS, sizeof(double));//change from 1 to DETECTOR_RESPONSE_BINS
	TempConst=(PROFILE_RESOLUTION)/(2*Geometry->delta_xz);
	OffsetR = ((Geometry->delta_xz/sqrt(3)) + Sinogram->delta_r/2)/DETECTOR_RESPONSE_BINS;
	OffsetT = ((Geometry->delta_xy/2) + Sinogram->delta_t/2)/DETECTOR_RESPONSE_BINS;

	for(k = 0 ; k < Sinogram->N_theta; k++)
	{
		for (i = 0; i < DETECTOR_RESPONSE_BINS; i++) //displacement along r
		{
			ProfileCenterR = i*OffsetR;
			rmin = ProfileCenterR - Geometry->delta_xz;
			for (j = 0 ; j < 1; j++)//displacement along t ;change to DETECTOR_RESPONSE_BINS later
			{
				ProfileCenterT = j*OffsetT;
				tmin = ProfileCenterT - Geometry->delta_xy/2;
				sum = 0;
				for (p=0; p < BEAM_RESOLUTION; p++)
				{
					r = r0 + p*StepSize;
					if(r < rmin)
						continue;

					ProfileIndex = (int32_t)floor((r-rmin)*TempConst);
					if(ProfileIndex < 0)
						ProfileIndex = 0;
					if(ProfileIndex >= PROFILE_RESOLUTION)
						ProfileIndex = PROFILE_RESOLUTION-1;

					//for (l =0; l < BEAM_RESOLUTION; l++)
					//{
					//t = t0 + l*StepSize;
					//if( t < tmin)
					//   continue;
					sum += (VoxelProfile[k][ProfileIndex] * BeamProfile[p]);//;*BeamProfile[l]);
					//}

				}
				H[j][k][i] = sum;

			}
		}
	}
	printf("Detector Done\n");
	fwrite(&H[0][0][0], sizeof(double),Sinogram->N_theta*DETECTOR_RESPONSE_BINS, Fp);
	fclose(Fp);
	return H;
}



#ifndef STORE_A_MATRIX
/*
void* CE_CalculateAMatrixColumnPartial(uint16_t row,uint16_t col,Sino* Sinogram,Geom* Geometry,double** VoxelProfile)
{
	int32_t i,j,k;
	double x,z,y;
	double r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
	double rmax,rmin;//stores the start and end points of the pixel profile on the detector
	double RTemp,TempConst,checksum = 0,Integral = 0;
	double LeftEndOfBeam;
	double MaximumSpacePerColumn;//we will use this to allocate space
	double AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults

	int32_t index_min,index_max;//stores the detector index in which the profile lies
	int32_t BaseIndex,FinalIndex,ProfileIndex=0;
	uint32_t count = 0;



	AMatrixCol* Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	AMatrixCol* Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end



	x = Geometry->x0 + ((double)col+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	z = Geometry->z0 + ((double)row+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner

	TempConst=(PROFILE_RESOLUTION)/(2*Geometry->delta_xz);



	AvgNumXElements = ceil(3*Geometry->delta_xz/Sinogram->delta_r);
	AvgNumYElements = ceil(3*Geometry->delta_xy/Sinogram->delta_t);
	MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*Sinogram->N_theta;

	Temp->values = (double*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(double));
	Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));

	//printf("%lf",Temp->values[10]);

#ifdef AREA_WEIGHTED
	for(i=0;i<Sinogram->N_theta;i++)
	{

		r = x*cosine[i] - z*sine[i];

		rmin = r - Geometry->delta_xz;
		rmax = r + Geometry->delta_xz;

		if(rmax < Sinogram->R0 || rmin > Sinogram->RMax)
			continue;



		index_min = floor(((rmin - Sinogram->R0)/Sinogram->delta_r));
		index_max = floor((rmax - Sinogram->R0)/Sinogram->delta_r);

		if(index_max >= Sinogram->N_r)
			index_max = Sinogram->N_r - 1;

		if(index_min < 0)
			index_min = 0;

		BaseIndex = i*Sinogram->N_r;

		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);


		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are
		//of same dimesions


		for(j = index_min;j <= index_max; j++)//Check
		{

			Integral = 0.0;

			//Accounting for Beam width
			RTemp = (Sinogram->R0 + (((double)j) + 0.5) *(Sinogram->delta_r));//the 0.5 is to get to the center of the detector

			LeftEndOfBeam = RTemp - (BEAM_WIDTH/2);//Consider pre storing the divide by 2

			//   if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
			//printf("%d %lf %lf\n",j, LeftEndOfBeam,rmin);


			for(k=0; k < BEAM_RESOLUTION; k++)
			{

				RTemp = LeftEndOfBeam + ((((double)k)*(BEAM_WIDTH))/BEAM_RESOLUTION);

				if (RTemp-rmin >= 0.0)
				{
					ProfileIndex = (uint32_t)floor((RTemp-rmin)*TempConst);//Finding the nearest neighbor profile to the beam
					//if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
					//printf("%d\n",ProfileIndex);
					if(ProfileIndex > PROFILE_RESOLUTION)
						ProfileIndex = PROFILE_RESOLUTION;
				}
				if(ProfileIndex < 0)
					ProfileIndex = 0;


				if(ProfileIndex >= 0 && ProfileIndex < PROFILE_RESOLUTION)
				{
#ifdef BEAM_CALCULATION
					Integral+=(BeamProfile[k]*VoxelProfile[i][ProfileIndex]);
#else
					Integral+=(VoxelProfile[i][ProfileIndex]/BEAM_RESOLUTION);
#endif
					//	if(FinalIndex >=176 && FinalIndex <= 178 && row ==0 && col == 0)
					//	 printf("Index %d %lf Voxel %lf I=%d\n",ProfileIndex,BeamProfile[k],VoxelProfile[2][274],i);
				}

			}
			if(Integral > 0.0)
			{
						FinalIndex = BaseIndex + (int32_t)j;
						Temp->values[count] = Integral;
						Temp->index[count] = FinalIndex;//can instead store a triple (row,col,slice) for the sinogram
						count++;


			}


		}



	}

#endif

	// printf("Final Space allocation for column %d %d\n",row,col);

	Ai->values=(double*)get_spc(count,sizeof(double));
	Ai->index=(uint32_t*)get_spc(count,sizeof(uint32_t));
	k=0;
	for(i = 0; i < count; i++)
	{
		if(Temp->values[i] > 0.0)
		{
			Ai->values[k]=Temp->values[i];
			checksum+=Ai->values[k];
			Ai->index[k++]=Temp->index[i];
		}

	}

	Ai->count=k;

	//printf("%d %d \n",Ai->count,count);

	//printf("(%d,%d) %lf \n",row,col,checksum);

	free(Temp->values);
	free(Temp->index);
	free(Temp);
	return Ai;

}
*/

void* CE_CalculateAMatrixColumnPartial(uint16_t row,uint16_t col, uint16_t slice, Sino* Sinogram,Geom* Geometry,double*** DetectorResponse)
{
	int32_t i,j,k,sliceidx;
	double x,z,y;
	double r;//this is used to find where does the ray passing through the voxel at certain angle hit the detector
	double t; //this is similar to r but along the y direction
	double tmin,tmax;
	double rmax,rmin;//stores the start and end points of the pixel profile on the detector
	double R_Center,TempConst,checksum = 0,Integral = 0,delta_r;
	double T_Center,delta_t;
	double MaximumSpacePerColumn;//we will use this to allocate space
	double AvgNumXElements,AvgNumYElements;//This is a measure of the expected amount of space per Amatrixcolumn. We will make a overestimate to avoid seg faults
	double ProfileThickness,stepsize;
	double OffsetR;
	double OffsetT;
	//interpolation variables
	double w1,w2,w3,w4,f1,f2,InterpolatedValue,ContributionAlongT;
	int32_t index_min,index_max,slice_index_min,slice_index_max,index_delta_r,index_delta_t;//stores the detector index in which the profile lies
	int32_t BaseIndex,FinalIndex,ProfileIndex=0;
	int32_t NumOfDisplacements=32;
	uint32_t count = 0;
  AMatrixCol* Ai = NULL;
  AMatrixCol* Temp = NULL;

	OffsetR = ((Geometry->delta_xz/sqrt(3)) + Sinogram->delta_r/2)/DETECTOR_RESPONSE_BINS;
	OffsetT = ((Geometry->delta_xz/2) + Sinogram->delta_t/2)/DETECTOR_RESPONSE_BINS;

	Ai = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	Temp = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));//This will assume we have a total of N_theta*N_x entries . We will freeuname -m this space at the end

	x = Geometry->x0 + ((double)col+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	z = Geometry->z0 + ((double)row+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	y = Geometry->y0 + ((double)slice + 0.5)*Geometry->delta_xy;

	TempConst=(PROFILE_RESOLUTION)/(2*Geometry->delta_xz);

	//alternately over estimate the maximum size require for a single AMatrix column
	AvgNumXElements = ceil(3*Geometry->delta_xz/Sinogram->delta_r);
	AvgNumYElements = ceil(3*Geometry->delta_xy/Sinogram->delta_t);
	MaximumSpacePerColumn = (AvgNumXElements * AvgNumYElements)*Sinogram->N_theta;

	Temp->values = (double*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(double));
	Temp->index  = (uint32_t*)get_spc((uint32_t)MaximumSpacePerColumn,sizeof(uint32_t));


#ifdef AREA_WEIGHTED
	for(i=0;i<Sinogram->N_theta;i++)
	{

		r = x*cosine[i] - z*sine[i];
		t = y;

		rmin = r - Geometry->delta_xz;
		rmax = r + Geometry->delta_xz;

		tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
		tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax ? t + Geometry->delta_xy/2 : Sinogram->TMax;

		if(rmax < Sinogram->R0 || rmin > Sinogram->RMax)
			continue;



		index_min = floor(((rmin - Sinogram->R0)/Sinogram->delta_r));
		index_max = floor((rmax - Sinogram->R0)/Sinogram->delta_r);

		if(index_max >= Sinogram->N_r)
			index_max = Sinogram->N_r - 1;

		if(index_min < 0)
			index_min = 0;

		slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_t);
		slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_t);

		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= Sinogram->N_t)
			slice_index_max = Sinogram->N_t -1;

		BaseIndex = i*Sinogram->N_r;//*Sinogram->N_t;

		for(j = index_min;j <= index_max; j++)//Check
		{



			//Accounting for Beam width
			R_Center = (Sinogram->R0 + (((double)j) + 0.5) *(Sinogram->delta_r));//the 0.5 is to get to the center of the detector

			//Find the difference between the center of detector and center of projection and compute the Index to look up into
			delta_r = fabs(r - R_Center);
			index_delta_r = floor((delta_r/OffsetR));


			if (index_delta_r >= 0 && index_delta_r < DETECTOR_RESPONSE_BINS)
			{

		//		for (sliceidx = slice_index_min; sliceidx <= slice_index_max; sliceidx++)
		//		{
					T_Center = (Sinogram->T0 + (((double)sliceidx) + 0.5) *(Sinogram->delta_t));
					delta_t = fabs(t - T_Center);
					index_delta_t = 0;//floor(delta_t/OffsetT);



					if (index_delta_t >= 0 && index_delta_t < DETECTOR_RESPONSE_BINS)
					{

						//Using index_delta_t,index_delta_t+1,index_delta_r and index_delta_r+1 do bilinear interpolation
						w1 = delta_r - index_delta_r*OffsetR;
						w2 = (index_delta_r+1)*OffsetR - delta_r;

						w3 = delta_t - index_delta_t*OffsetT;
						w4 = (index_delta_r+1)*OffsetT - delta_t;



						f1 = (w2/OffsetR)*DetectorResponse[index_delta_t][i][index_delta_r] + (w1/OffsetR)*DetectorResponse[index_delta_t][i][index_delta_r+1 < DETECTOR_RESPONSE_BINS ? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];
						//	f2 = (w2/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS ?index_delta_t+1 : DETECTOR_RESPONSE_BINS-1][i][index_delta_r] + (w1/OffsetR)*DetectorResponse[index_delta_t+1 < DETECTOR_RESPONSE_BINS? index_delta_t+1:DETECTOR_RESPONSE_BINS][i][index_delta_r+1 < DETECTOR_RESPONSE_BINS? index_delta_r+1:DETECTOR_RESPONSE_BINS-1];

						if(sliceidx == slice_index_min)
							ContributionAlongT = (sliceidx + 1)*Sinogram->delta_t - tmin;
						else if(sliceidx == slice_index_max)
							ContributionAlongT = tmax - (sliceidx)*Sinogram->delta_t;
						else {
							ContributionAlongT = Sinogram->delta_t;
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


	Ai->values=(double*)get_spc(count,sizeof(double));
	Ai->index=(uint32_t*)get_spc(count,sizeof(uint32_t));
	k=0;
	for(i = 0; i < count; i++)
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

#endif


