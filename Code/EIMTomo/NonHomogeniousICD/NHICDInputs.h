#ifndef NCICD_COMPUTATIONINPUTS_H_
#define NCICD_COMPUTATIONINPUTS_H_

#include <stdio.h> //For all other declarations int,FILE etc

#include "TomoEngine/TomoEngine.h"


#ifdef __cplusplus
extern "C" {
#endif

	/* Axes conventions:

	      . Y
	     .
	    .
	   .
	  .
	 .
	 ---------------------> X
	 |
	 |
	 |
	 |
	 |
	 |
	 V
	 Z
	 */

	struct _sinogram
	{
		uint16_t N_r;//Number of measurements in x direction
		uint16_t N_theta;//Number of angles
		uint16_t N_t;//Number of measurements in y direction
		double delta_r;//Distance between successive measurements along x
		double delta_t;//Distance between successive measurements along y
		double ***counts;//The measured images should be stored in this once read from the input file. It will be a Ny X (Nz X Nx)
		double *angles;//Holds the angles through which the object is tilted
		double R0,RMax;
		double T0,TMax;

	};

	typedef struct _sinogram Sinogram;

	struct _geometry
	{
		//User Input
		double LengthZ;//This is the sample thickness
		double delta_xz;//Voxel size in the x-z plane (assuming square shaped voxels in the x-z plane)
		double delta_xy;//Voxel size in the x-y plane
		double ***Object;//Holds the volume to be reconstructed

		//Computed From User Input
		double LengthX;//sinogram.N_x * delta_r;
		double LengthY;//sinogram.N_y * delta_t
		uint16_t N_x;//Number of voxels in x direction
		uint16_t N_z;//Number of voxels in z direction
		uint16_t N_y;//Number of measurements in y direction
		//Coordinates of the left corner of the x-z object
		double x0;// -LengthX/2
		double z0;// -LengthZ/2
		double y0;//-LengthY/2

	};

	typedef struct _geometry Geometry;

	struct _command_line_inputs
	{
		char* ParamFile;
		char* SinoFile;
		char* InitialRecon;
		char* OutputFile;
		//This is read from the paramter file
		int16_t NumIter;
		double SigmaX;
		double p;
	};
	typedef struct _command_line_inputs TomoInputs;


	void readParameterFile(FILE* ,TomoInputs* ,Sinogram* ,Geometry*);
	void initializeSinoParameters(Sinogram *,TomoInputs*);
	void initializeGeomParameters(Sinogram* ,Geometry* ,TomoInputs*);


#ifdef __cplusplus
}
#endif

#endif //ComputationInputs
