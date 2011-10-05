/*
 *  ComputationEngineOptimized.h
 *  HAADFSTEM_Xcode
 *
 *  Created by Singanallur Venkatakrishnan on 6/29/11.
 *  Copyright 2011 Purdue University. All rights reserved.
 *
 */

#ifndef MOTION_CORRECTION_COMPUTATIONENGINE_H_
#define MOTION_CORRECTION_COMPUTATIONENGINE_H_

#include "EIMTomo/EIMTomo.h"
#include "EIMTomo/common/EMMath.h"
#include "MotionCorrectionInputs.h"




	//#define DEBUG
#define PROFILE_RESOLUTION 1536
#define PI M_PI
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



#ifdef __cplusplus
  extern "C" {
#endif



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

