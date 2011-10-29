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
	
#include "ScaleOffsetCorrectionInputs.h"
	
	
	//#define DEBUG	,
#define PROFILE_RESOLUTION 1536
//#define PI 4*atan(1)//3.14159265
	//Beam Parameters - This is set to some number <<< Sinogram->delta_r.
//#define BEAM_WIDTH 0.050000 
#define BEAM_RESOLUTION 512
#define AREA_WEIGHTED
#define ROI //Region Of Interest
#define STOPPING_THRESHOLD 0.009
//#define	SURROGATE_FUNCTION
	//#define DISTANCE_DRIVEN
	//#define CORRECTION
//#define STORE_A_MATRIX
//	#define WRITE_INTERMEDIATE_RESULTS
#define COST_CALCULATE
//#define BEAM_CALCULATION
#define DETECTOR_RESPONSE_BINS 64
#define JOINT_ESTIMATION
//#define NOISE_MODEL
//#define	CIRCULAR_BOUNDARY_CONDITION
	
	//Structure to store a single column(A_i) of the A-matrix 
	typedef struct 
	{
		DATA_TYPE* values;//Store the non zero entries
		uint32_t count;//The number of non zero values present in the column
		uint32_t *index;//This maps each value to its location in the column. The entries in this can vary from 0 to Sinogram.N_x Sinogram.N_theta-1
	}AMatrixCol;
	
	typedef struct
	{
		DATA_TYPE *I_0;//Scale
		DATA_TYPE *mu;//Offset
	}ScaleOffsetParams;
	
	DATA_TYPE*** ForwardProject(Sino*,Geom*,DATA_TYPE***,DATA_TYPE***);
	void CE_CalculateSinCos(Sino*);
	void CE_InitializeBeamProfile(Sino*);
	void CE_MinMax(DATA_TYPE*,DATA_TYPE*);
	void* CE_CalculateVoxelProfile(Sino*,Geom*);
	int CE_MAPICDReconstruct(Sino*,Geom*,CommandLineInputs*);
	void* CE_CalculateAMatrixColumn(uint16_t ,uint16_t , uint16_t ,Sino* ,Geom* ,DATA_TYPE**);
	double CE_DerivOfCostFunc(double);
	DATA_TYPE CE_ComputeCost(DATA_TYPE***  ,DATA_TYPE***  ,Sino* ,Geom* );
	void* CE_CalculateAMatrixColumnPartial(uint16_t,uint16_t,uint16_t,Sino*,Geom*,DATA_TYPE***);
//	void* CE_CalculateAMatrixColumnPartial(uint16_t,uint16_t,Sino*,Geom*,DATA_TYPE**);
	void* CE_DetectorResponse(uint16_t ,uint16_t ,Sino* ,Geom* ,DATA_TYPE**);
	double  solve(
				  double (*f)(), /* pointer to function to be solved */
				  double a,      /* minimum value of solution */
				  double b,      /* maximum value of solution */
				  double err,    /* accuarcy of solution */
				  int *code      /* error code */
				  );
	double CE_SurrogateFunctionBasedMin();
	double CE_ConstraintEquation(DATA_TYPE);
	DATA_TYPE* CE_RootsOfQuadraticFunction(DATA_TYPE,DATA_TYPE,DATA_TYPE);
	
#ifdef __cplusplus
}
#endif

#endif //CompEngine

