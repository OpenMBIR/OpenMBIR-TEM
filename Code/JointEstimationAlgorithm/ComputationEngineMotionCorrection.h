/*
 *  ComputationEngineOptimized.h
 *  HAADFSTEM_Xcode
 *
 *  Created by Singanallur Venkatakrishnan on 6/29/11.
 *  Copyright 2011 Purdue University. All rights reserved.
 *
 */

#ifndef COMPUTATIONENGINE_H_
#define COMPUTATIONENGINE_H_
#ifdef __cplusplus
extern "C" {
#endif
	
#include "ComputationInputsMotionCorrection.h"
	
	
	//#define DEBUG	
#define PROFILE_RESOLUTION 1536
#define PI 4*atan(1)//3.14159265
	//Beam Parameters - This is set to some number <<< Sinogram->delta_r.
	//#define BEAM_WIDTH 0.050000 
#define BEAM_RESOLUTION 512
#define AREA_WEIGHTED
	//#define ROI //Region Of Interest
	//#define DISTANCE_DRIVEN
	//#define CORRECTION
	//#define STORE_A_MATRIX
	//	#define WRITE_INTERMEDIATE_RESULTS
	//#define COST_CALCULATE
#define BEAM_CALCULATION
#define DETECTOR_RESPONSE_BINS 64

#ifdef CORRECTION
	double *NORMALIZATION_FACTOR; 
#endif
	//Structure to store a single column(A_i) of the A-matrix 
	typedef struct 
	{
		double* values;//Store the non zero entries
		uint32_t count;//The number of non zero values present in the column
		uint32_t *index;//This maps each value to its location in the column. The entries in this can vary from 0 to Sinogram.N_x Sinogram.N_theta-1
	}AMatrixCol;
	
	
	//Markov Random Field Prior parameters - Globals -:( 
	double FILTER[3][3][3]={{{0.0302,0.0370,0.0302},{0.0370,0.0523,0.0370},{0.0302,0.0370,0.0302}},
		{{0.0370,0.0523,0.0370},{0.0523,0.0,0.0523},{0.0370,0.0523,  0.0370}},
		{{0.0302,0.0370,0.0302},{0.0370,0.0523,0.0370},{0.0302,0.0370,0.0302}}};
	double THETA1,THETA2,NEIGHBORHOOD[3][3][3],V;
	double MRF_P ;//= 1.1;
	double SIGMA_X_P;
	
	//Global variables 
	double *cosine,*sine;//used to store cosine and sine of all angles through which sample is tilted
	double *BeamProfile;//used to store the shape of the e-beam
    double BEAM_WIDTH;
	double OffsetR;
	double OffsetT;
	
	
	
	void CE_CalculateSinCos(Sino*);
	void CE_InitializeBeamProfile(Sino*);
	void CE_MinMax(double*,double*);
	void* CE_CalculateVoxelProfile(Sino*,Geom*);
	int CE_MAPICDReconstruct(Sino*,Geom*,CommandLineInputs*);
	void* CE_CalculateAMatrixColumn(uint16_t ,uint16_t , uint16_t ,Sino* ,Geom* ,double**);
	double CE_DerivOfCostFunc(double);
	double CE_ComputeCost(double***  ,double***  ,Sino* ,Geom* );
	void* CE_CalculateAMatrixColumnPartial(uint16_t,uint16_t,uint16_t,Sino*,Geom*,double***);
	//	void* CE_CalculateAMatrixColumnPartial(uint16_t,uint16_t,Sino*,Geom*,double**);
	void* CE_DetectorResponse(uint16_t ,uint16_t ,Sino* ,Geom* ,double**);
	double  solve(
				  double (*f)(), /* pointer to function to be solved */
				  double a,      /* minimum value of solution */
				  double b,      /* maximum value of solution */
				  double err,    /* accuarcy of solution */
				  int *code      /* error code */
				  );
	
#ifdef __cplusplus
}
#endif

#endif //CompEngine

