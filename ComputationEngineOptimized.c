/*
 *  ComputationEngineOptimized.c
 *  HAADFSTEM_Xcode
 *
 *  Created by Singanallur Venkatakrishnan on 6/29/11.
 *  Copyright 2011 Purdue University. All rights reserved.
 *
 */

#include "ComputationEngineOptimized.h"

#include "allocate.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//#define DEBUG


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
	int16_t Iter;
	int16_t i,j,k,r,row,col,slice,RowIndex,ColIndex,SliceIndex,Idx;
	int16_t NumOfXPixels;
	int32_t q,p;
	double checksum = 0,temp;
	double *cost;
	double** VoxelProfile;
	double ***Y_Est;//Estimted Sinogram 
	double ***ErrorSino;//Error Sinogram
	double ***Weight;//This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation
#ifdef DEBUG
	FILE *Fp = fopen("ReconstructedSino.bin","w");//Reconstructed Sinogram from initial est
#endif
	FILE* Fp2;//Cost function 
	FILE *Fp3;//File to store intermediate outputs of reconstruction
	
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
	
	Y_Est=(double ***)get_3D(Sinogram->N_y,Sinogram->N_theta,Sinogram->N_x,sizeof(double));
	ErrorSino=(double ***)get_3D(Sinogram->N_y,Sinogram->N_theta,Sinogram->N_x,sizeof(double));
	Weight=(double ***)get_3D(Sinogram->N_y,Sinogram->N_theta,Sinogram->N_x,sizeof(double));
	
	//calculate the trapezoidal voxel profile for each angle.Also the angles in the Sinogram Structure are converted to radians
	VoxelProfile = CE_CalculateVoxelProfile(Sinogram,Geometry); //Verified with ML
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
			if (((x*x)/(EllipseA*EllipseA)) + ((z*z)/(EllipseB*EllipseB)) < 0.95)  
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
	NORMALIZATION_FACTOR = (double*)get_spc(Sinogram->N_x*Sinogram->N_theta,sizeof(double));
	
	for(i = 0;i < Sinogram->N_x*Sinogram->N_theta;i++)
		NORMALIZATION_FACTOR[i] = 1.0;
	
#endif 
	
	//Assign the scaling to a global variable
#ifdef BEAM_CALCULATION
	BEAM_WIDTH = (1)*Sinogram->delta_x;
#else
	BEAM_WIDTH = Sinogram->delta_x;
#endif
	
	MRF_P = CmdInputs->p;
	SIGMA_X_P = pow(CmdInputs->SigmaX,MRF_P);
	
	for(i=0;i<Sinogram->N_theta;i++)
		for(j=0;j< PROFILE_RESOLUTION;j++)
			checksum+=VoxelProfile[i][j];
    printf("CHK SUM%lf\n",checksum);
	
    checksum=0;
    
	
	//Initialize the e-beam 
	CE_InitializeBeamProfile(Sinogram); //verified with ML
	
	//calculate sine and cosine of all angles and store in the global arrays sine and cosine
	
	CE_CalculateSinCos(Sinogram);
	
	
	//Allocate space for storing columns the A-matrix; an array of pointers to columns 
	//AMatrixCol** AMatrix=(AMatrixCol **)get_spc(Geometry->N_x*Geometry->N_z,sizeof(AMatrixCol*));
	
#ifdef STORE_A_MATRIX
	
	AMatrixCol**** AMatrix = (AMatrixCol ****)multialloc(sizeof(AMatrixCol*),3,Geometry->N_y,Geometry->N_z,Geometry->N_x);
#else
	double y;
	double t,tmin,tmax,ProfileThickness; 
	int16_t slice_index_min,slice_index_max;
	AMatrixCol*** TempCol = (AMatrixCol***)multialloc(sizeof(AMatrixCol*),2,Geometry->N_z,Geometry->N_x);
#endif
	
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
    TempCol[j][k] = CE_CalculateAMatrixColumnPartial(j,k,Sinogram,Geometry,VoxelProfile);
    }
#endif
	
	printf("Geometry-Z %d\n",Geometry->N_z);
	
	//Forward Project Geometry->Object one slice at a time and compute the  Sinogram for each slice
	//is Y_Est initailized to zero?
	for(i = 0; i < Sinogram->N_y; i++)
		for(j = 0; j < Sinogram->N_theta; j++)
			for(k = 0;k < Sinogram->N_x; k++)
				Y_Est[i][j][k]=0.0;
	
	
	
	
	for (i = 0;i < Geometry->N_y; i++)//slice index
	{
		//  p=0;//This is used to access the A-Matrix column correspoding the (i,j,k)the voxel location
		for(j = 0;j < Geometry->N_z; j++)
			for(k = 0;k < Geometry->N_x; k++)
			{
#ifdef STORE_A_MATRIX  //A matrix call
				for(q = 0;q < AMatrix[i][j][k]->count; q++)
					//for(q = 0; q < AMatrix[j][k]->count; q++)
				{
					//Convert AMatrix->index[p] to a (slice,row,col)
					row = (int16_t)floor(AMatrix[i][j][k]->index[q]/(Sinogram->N_x*Sinogram->N_y));
					slice =(int16_t) floor((AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))/(Sinogram->N_x));
					col =  (int16_t)(AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))%(Sinogram->N_x);
					//row = (uint32_t)((AMatrix[j][k]->index[q])/Sinogram->N_x);
					//col = (AMatrix[j][k]->index[q])%Sinogram->N_x;
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
					
		           slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_y);
		           slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_y);
		
		           if(slice_index_min < 0)
			       slice_index_min = 0;
		           if(slice_index_max >= Sinogram->N_y)
			       slice_index_max = Sinogram->N_y-1; 

				for(q = 0;q < TempCol[j][k]->count; q++)
				{
				    //calculating the footprint of the voxel in the t-direction
				    			    					
					RowIndex = floor(TempCol[j][k]->index[q]/(Sinogram->N_x));
					ColIndex =  (TempCol[j][k]->index[q]%(Sinogram->N_x));
					for(SliceIndex = slice_index_min ; SliceIndex <= slice_index_max; SliceIndex++)
					{
					if(SliceIndex == slice_index_min)
						ProfileThickness = (((SliceIndex+1)*Sinogram->delta_y + Sinogram->T0) - tmin);//Sinogram->delta_y; //Will be < Sinogram->delta_y
					else 
					{
						if (SliceIndex == slice_index_max) 
						{
							ProfileThickness = (tmax - ((SliceIndex)*Sinogram->delta_y + Sinogram->T0));//Sinogram->delta_y;//Will be < Sinogram->delta_y	
						}
						else 
						{
							ProfileThickness = Sinogram->delta_y;
						}
						
					}
					Y_Est[SliceIndex][RowIndex][ColIndex] += (TempCol[j][k]->values[q] *ProfileThickness* Geometry->Object[i][j][k]);
					}
					
				}
				
#endif
			} 
		
	}
	
	
	//Calculate Error Sinogram - Can this be combined with previous loop?
	//Also compute weights of the diagonal covariance matrix
	
	for (i=0;i<Sinogram->N_y;i++)//slice index
	{
		checksum=0;
		for(j=0;j<Sinogram->N_theta;j++)
			for(k=0;k<Sinogram->N_x;k++)
			{
				
				// checksum+=Y_Est[i][j][k];
				ErrorSino[i][j][k] = Sinogram->counts[i][j][k] - Y_Est[i][j][k];
				if(Sinogram->counts[i][j][k] != 0)
					Weight[i][j][k]=1.0/Sinogram->counts[i][j][k];
				else
					Weight[i][j][k]=0;
#ifdef DEBUG			
				//writing the error sinogram
				fwrite(&ErrorSino[i][j][k],sizeof(double),1,Fp);
#endif	
				checksum+=Weight[i][j][k];
			}
		printf("Check sum = %lf\n",checksum);  
		
	}
	
	free(Y_Est);
#ifdef DEBUG
	fclose(Fp);
#endif
	
	
	
	
#ifdef COST_CALCULATE
	Fp2 = fopen("CostFunc.bin","w");
	
	printf("Cost Function Calcluation \n");
	/*********************Initial Cost Calculation***************************************************/
	temp=0;
	cost[0]=0;
	for (i = 0; i < Sinogram->N_y; i++)
		for (j = 0; j < Sinogram->N_theta; j++)
			for( k = 0; k < Sinogram->N_x; k++)
				cost[0] += (ErrorSino[i][j][k] * ErrorSino[i][j][k] * Weight[i][j][k]);
	cost[0]/=2;//accounts for the half in from: (1/2)||Y-AX||^2 + ...
	
	//Each neighboring pair should be counted only once
	
	for (i = 0; i < Geometry->N_y; i++)
		for (j = 0; j < Geometry->N_z; j++)
			for(k = 0; k < Geometry->N_x; k++)
			{
				
				if(k+1 <  Geometry->N_x)
					temp += FILTER[1][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j][k+1]),MRF_P);
				
				
				if(j+1 < Geometry->N_z)
				{
					if(k-1 >= 0)
						temp += FILTER[1][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k-1]),MRF_P);
					
					
					temp += FILTER[1][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k]),MRF_P);
					
					
					if(k+1 < Geometry->N_x)
						temp += FILTER[1][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k+1]),MRF_P);
					
				}
				
				if(i+1 < Geometry->N_y)
				{
					
					if(j-1 >= 0)
						temp += FILTER[2][0][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k]),MRF_P);
					
					temp += FILTER[2][1][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k]),MRF_P);
					
					if(j+1 < Geometry->N_z)
						temp += FILTER[2][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k]),MRF_P);
					
					
					if(j-1 >= 0)
					{
						if(k-1 >= 0)
							temp += FILTER[2][0][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k-1]),MRF_P);
						
						if(k+1 < Geometry->N_x)
							temp += FILTER[2][0][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k+1]),MRF_P);
						
					}
					
					if(k-1 >= 0)
						temp += FILTER[2][1][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k-1]),MRF_P);
					
					if(j+1 < Geometry->N_z)
					{
						if(k-1 >= 0)
							temp += FILTER[2][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k-1]),MRF_P);
						
						if(k+1 < Geometry->N_x)
							temp+= FILTER[2][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k+1]),MRF_P);
					}
					
					if(k+1 < Geometry->N_x)
						temp+= FILTER[2][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k+1]),MRF_P);
				}
			}
	
	cost[0] += (temp/(MRF_P*SIGMA_X_P));
	
	printf("%lf\n",cost[0]);
	
	fwrite(&cost[0],sizeof(double),1,Fp2);
	/*******************************************************************************/
#endif //Cost calculation endif
	
	
	//Loop through every voxel updating it by solving a cost function
	
	for(Iter = 0;Iter < CmdInputs->NumIter;Iter++)
	{
		
		
		//printf("Iter %d\n",Iter);
		for (i = 0; i < Geometry->N_y; i++)//slice index
		{
			//Idx=0;
			//printf("Slice %d\n",i);
			for(j = 0; j < Geometry->N_z; j++)//Row index
				for(k = 0; k < Geometry->N_x ; k++)//Column index
				{
#ifdef ROI
					if (Mask[j][k] == 1)
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
										if(j+q >= 0 && j+q < Geometry->N_z)
											if(k+r >= 0 && k+r < Geometry->N_x)
											{
												NEIGHBORHOOD[p+1][q+1][r+1] = Geometry->Object[p+i][q+j][r+k];
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
						
						V = Geometry->Object[i][j][k];//Store the present value of the voxel
						THETA1 = 0.0;
						THETA2 = 0.0;
						
#ifdef STORE_A_MATRIX
						
						for (q = 0;q < AMatrix[i][j][k]->count; q++)
						{
							
							
							RowIndex = floor(AMatrix[i][j][k]->index[q]/(Sinogram->N_x*Sinogram->N_y));
							ColIndex = (AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))%(Sinogram->N_x);
							SliceIndex = floor((AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))/(Sinogram->N_x)); 
							
							THETA2 += (AMatrix[i][j][k]->values[q])*(AMatrix[i][j][k]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
							THETA1 += ErrorSino[SliceIndex][RowIndex][ColIndex]*(AMatrix[i][j][k]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
						}
#else
						
						y = ((double)i+0.5)*Geometry->delta_xy + Geometry->y0;
						t = y;
						tmin = (t - Geometry->delta_xy/2) > Sinogram->T0 ? t-Geometry->delta_xy/2 : Sinogram->T0;
						tmax = (t + Geometry->delta_xy/2) <= Sinogram->TMax? t + Geometry->delta_xy/2 : Sinogram->TMax;
						
						slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_y);
						slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_y);
						
						if(slice_index_min < 0)
							slice_index_min = 0;
						if(slice_index_max >= Sinogram->N_y)
							slice_index_max = Sinogram->N_y-1;
						
						//TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
						for (q = 0;q < TempCol[j][k]->count; q++)
						{
							
							RowIndex = floor(TempCol[j][k]->index[q]/(Sinogram->N_x));
							ColIndex =  (TempCol[j][k]->index[q]%(Sinogram->N_x));
							 for(SliceIndex = slice_index_min ; SliceIndex <= slice_index_max; SliceIndex++)
							 {
								 if(SliceIndex == slice_index_min)
									 ProfileThickness = (((SliceIndex+1)*Sinogram->delta_y + Sinogram->T0) - tmin);//Sinogram->delta_y; //Will be < Sinogram->delta_y
								 else 
								 {
									 if (SliceIndex == slice_index_max) 
									 {
										 ProfileThickness = (tmax - ((SliceIndex)*Sinogram->delta_y + Sinogram->T0));//Sinogram->delta_y;//Will be < Sinogram->delta_y	
									 }
									 else 
									 {
										 ProfileThickness = Sinogram->delta_y;
									 }
									 
								 }
								 
								 
								 THETA2 += (ProfileThickness*ProfileThickness)*(TempCol[j][k]->values[q])*(TempCol[j][k]->values[q])*Weight[SliceIndex][RowIndex][ColIndex];
								 THETA1 += ErrorSino[SliceIndex][RowIndex][ColIndex]*(TempCol[j][k]->values[q])*(ProfileThickness)*Weight[SliceIndex][RowIndex][ColIndex];
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
						Geometry->Object[i][j][k] = UpdatedVoxelValue;
						//Update the ErrorSinogram
						
#ifdef STORE_A_MATRIX
						for (q = 0;q < AMatrix[i][j][k]->count; q++)
						{
							
							RowIndex = floor(AMatrix[i][j][k]->index[q]/(Sinogram->N_x*Sinogram->N_y));
							ColIndex = (AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))%(Sinogram->N_x);
							SliceIndex = floor((AMatrix[i][j][k]->index[q]%(Sinogram->N_x*Sinogram->N_y))/(Sinogram->N_x)); 
							
							ErrorSino[SliceIndex][RowIndex][ColIndex] -= (AMatrix[i][j][k]->values[q] * (Geometry->Object[i][j][k] - V));
						}
#else
							       
					for (q = 0;q < TempCol[j][k]->count; q++)
					{
							
					RowIndex = floor(TempCol[j][k]->index[q]/(Sinogram->N_x));
					ColIndex =  (TempCol[j][k]->index[q]%(Sinogram->N_x));
					for(SliceIndex = slice_index_min ; SliceIndex <= slice_index_max; SliceIndex++)
					{
					if(SliceIndex == slice_index_min)
						ProfileThickness = (((SliceIndex+1)*Sinogram->delta_y + Sinogram->T0) - tmin);//(Sinogram->delta_y); //Will be < Sinogram->delta_y
					else 
					{
						if (SliceIndex == slice_index_max) 
						{
							ProfileThickness = (tmax - ((SliceIndex)*Sinogram->delta_y + Sinogram->T0));//(Sinogram->delta_y);//Will be < Sinogram->delta_y	
						}
						else 
						{
							ProfileThickness = (Sinogram->delta_y);
						}
						
					}
											
							ErrorSino[SliceIndex][RowIndex][ColIndex] -= (TempCol[j][k]->values[q] *ProfileThickness* (Geometry->Object[i][j][k] - V));
					}
					}
#endif
						Idx++;
					}
					
#ifdef ROI
				} //closing the if for region of interest only updates 
#endif
		}
		
		if (CE_Cancel == 1)
		{
			return err;
		}
		
#ifdef COST_CALCULATE
		/*********************Cost Calculation***************************************************/
		temp=0;
		cost[Iter+1]=0;
		for (i = 0; i < Sinogram->N_y; i++)
			for (j = 0; j < Sinogram->N_theta; j++)
				for( k = 0; k < Sinogram->N_x; k++)
					cost[Iter+1] += (ErrorSino[i][j][k] * ErrorSino[i][j][k] * Weight[i][j][k]);
		cost[Iter+1]/=2;//Accounting for (1/2) in the cost function
		//Each neighboring pair should be counted only once
		
		for (i = 0; i < Geometry->N_y; i++)
			for (j = 0; j < Geometry->N_z; j++)
				for(k = 0; k < Geometry->N_x; k++)
				{
					
					if(k+1 <  Geometry->N_x)
						temp += FILTER[1][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j][k+1]),MRF_P);
					
					
					if(j+1 < Geometry->N_z)
					{
						if(k-1 >= 0)
							temp += FILTER[1][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k-1]),MRF_P);
						
						
						temp += FILTER[1][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k]),MRF_P);
						
						
						if(k+1 < Geometry->N_x)
							temp += FILTER[1][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i][j+1][k+1]),MRF_P);
						
					}
					
					if(i+1 < Geometry->N_y)
					{
						
						if(j-1 >= 0)
							temp += FILTER[2][0][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k]),MRF_P);
						
						temp += FILTER[2][1][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k]),MRF_P);
						
						if(j+1 < Geometry->N_z)
							temp += FILTER[2][2][1]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k]),MRF_P);
						
						
						if(j-1 >= 0)
						{
							if(k-1 >= 0)
								temp += FILTER[2][0][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k-1]),MRF_P);
							
							if(k+1 < Geometry->N_x)
								temp += FILTER[2][0][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j-1][k+1]),MRF_P);
							
						}
						
						if(k-1 >= 0)
							temp += FILTER[2][1][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k-1]),MRF_P);
						
						if(j+1 < Geometry->N_z)
						{
							if(k-1 >= 0)
								temp += FILTER[2][2][0]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k-1]),MRF_P);
							
							if(k+1 < Geometry->N_x)
								temp+= FILTER[2][2][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j+1][k+1]),MRF_P);
						}
						
						if(k+1 < Geometry->N_x)
							temp+= FILTER[2][1][2]*pow(fabs(Geometry->Object[i][j][k]-Geometry->Object[i+1][j][k+1]),MRF_P);
					}
				}
		
		cost[Iter] += (temp/(MRF_P*SIGMA_X_P));
		
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
	
	//Temp->index = (uint32_t*)get_spc(Sinogram->N_x*Sinogram->N_theta,sizeof(uint32_t));
	//Temp->values = (double*)multialloc(sizeof(double),1,Sinogram->N_x*Sinogram->N_theta);//makes the values =0
	
	x = Geometry->x0 + ((double)col+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	z = Geometry->z0 + ((double)row+0.5)*Geometry->delta_xz;//0.5 is for center of voxel. x_0 is the left corner
	y = Geometry->y0 + ((double)slice + 0.5)*Geometry->delta_xy;
	
	TempConst=(PROFILE_RESOLUTION)/(2*Geometry->delta_xz);
	
	
	//	Temp->values = (double*)calloc(Sinogram->N_y*Sinogram->N_x*Sinogram->N_theta,sizeof(double));//(double*)get_spc(Sinogram->N_x*Sinogram->N_theta,sizeof(double));//makes the values =0
	
	//alternately over estimate the maximum size require for a single AMatrix column
	AvgNumXElements = ceil(3*Geometry->delta_xz/Sinogram->delta_x);
	AvgNumYElements = ceil(3*Geometry->delta_xy/Sinogram->delta_y);
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
		
		
		
		index_min = floor(((rmin - Sinogram->R0)/Sinogram->delta_x));
		index_max = floor((rmax - Sinogram->R0)/Sinogram->delta_x);
		
		if(index_max >= Sinogram->N_x)
			index_max = Sinogram->N_x - 1;
		
		if(index_min < 0)
			index_min = 0;
		
		slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_y);
		slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_y);
		
		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= Sinogram->N_y)
			slice_index_max = Sinogram->N_y -1; 
		
		BaseIndex = i*Sinogram->N_x*Sinogram->N_y;
		
		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);
		
		
		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are 
		//of same dimesions
		
		
		for(j = index_min;j <= index_max; j++)//Check 
		{
			
			Integral = 0.0;
			
			//Accounting for Beam width
			RTemp = (Sinogram->R0 + (((double)j) + 0.5) *(Sinogram->delta_x));//the 0.5 is to get to the center of the detector
			
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
						ProfileThickness = (((sliceidx+1)*Sinogram->delta_y + Sinogram->T0) - tmin);//Sinogram->delta_y; //Will be < Sinogram->delta_y
					else 
					{
						if (sliceidx == slice_index_max) 
						{
							ProfileThickness = (tmax - ((sliceidx)*Sinogram->delta_y + Sinogram->T0));//Sinogram->delta_y;//Will be < Sinogram->delta_y	
						}
						else 
						{
							ProfileThickness = Sinogram->delta_y;//Sinogram->delta_y;
						}
						
					}
					if(ProfileThickness > 0)
					{
						
						FinalIndex = BaseIndex + (int32_t)j + (int32_t)sliceidx * Sinogram->N_x;
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
		
		
		if(rmax < Sinogram->R0 || rmin > Sinogram->R0 + Sinogram->N_x*Sinogram->delta_x)
			continue;
		
		index_min = floor(((rmin-Sinogram->R0)/Sinogram->delta_x));
		index_max = floor((rmax-Sinogram->R0)/Sinogram->delta_x);
		
		slice_index_min = floor((tmin - Sinogram->T0)/Sinogram->delta_y);
		slice_index_max = floor((tmax - Sinogram->T0)/Sinogram->delta_y);
		
		if(slice_index_min < 0)
			slice_index_min = 0;
		if(slice_index_max >= Sinogram->N_y)
			slice_index_max = Sinogram->N_y -1; 
		
		
		if(index_max >= Sinogram->N_x)
			index_max = Sinogram->N_x-1;
		
		if(index_min < 0)
			index_min = 0;
		
		BaseIndex = i*Sinogram->N_x*Sinogram->N_y;
		
		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);
		
		
		
		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are 
		//of same dimesions
		
		for(j = index_min;j <= index_max; j++)//Check 
		{
			d1  = ((double)j)*Sinogram->delta_x + Sinogram->R0;
			d2 =  d1 + Sinogram->delta_x;
			
			if(rmax < d1)
			{
				Integral = 0;
			}
			else
			{
				if(rmin > d1 && rmin < d2 && rmax > d2)
				{
					Integral = (d2 - rmin)/Sinogram->delta_x;
				}
				else
				{
					if(rmin >= d1 && rmin <= d2 && rmax >= d1 && rmax <= d2)
						Integral= (rmax - rmin)/Sinogram->delta_x;
					else
						if(rmax > d1 && rmax < d2 && rmin < d1)
							Integral = (rmax - d1)/Sinogram->delta_x;
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
						ProfileThickness = (((sliceidx+1)*Sinogram->delta_y + Sinogram->T0) - tmin)*Sinogram->delta_y;
					else {
						if (sliceidx == slice_index_max) 
						{
							ProfileThickness = (tmax - ((sliceidx)*Sinogram->delta_y + Sinogram->T0))*Sinogram->delta_y;	
						}
						else {
							ProfileThickness = Sinogram->delta_y;
						}
						
					}
					if (ProfileThickness > 0)
					{
						FinalIndex = BaseIndex + (uint32_t)j + (int32_t)sliceidx * Sinogram->N_x;
						
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
		BeamProfile[i] = (1.0/(BEAM_WIDTH)) * ( 1 + cos ((PI/W)*fabs(-W + i*(BEAM_WIDTH/BEAM_RESOLUTION))));
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
	
	for (i = 0; i < Sinogram->N_y; i++)
		for (j = 0; j < Sinogram->N_theta; j++)
			for( k = 0; k < Sinogram->N_x; k++)
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
#ifndef STORE_A_MATRIX
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
	
	
	
	AvgNumXElements = ceil(3*Geometry->delta_xz/Sinogram->delta_x);
	AvgNumYElements = ceil(3*Geometry->delta_xy/Sinogram->delta_y);
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
		
		
		
		index_min = floor(((rmin - Sinogram->R0)/Sinogram->delta_x));
		index_max = floor((rmax - Sinogram->R0)/Sinogram->delta_x);
		
		if(index_max >= Sinogram->N_x)
			index_max = Sinogram->N_x - 1;
		
		if(index_min < 0)
			index_min = 0;
	
		BaseIndex = i*Sinogram->N_x;
		
		// if(row == 63 && col == 63)
		//    printf("%d %d\n",index_min,index_max);
		
		
		//Do calculations only if there is a non zero contribution to the slice. This is particulary useful when the detector and geometry are 
		//of same dimesions
		
		
		for(j = index_min;j <= index_max; j++)//Check 
		{
			
			Integral = 0.0;
			
			//Accounting for Beam width
			RTemp = (Sinogram->R0 + (((double)j) + 0.5) *(Sinogram->delta_x));//the 0.5 is to get to the center of the detector
			
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
#endif