#ifndef COMPUTATIONINPUTS_H_
#define COMPUTATIONINPUTS_H_

#include <stdio.h> //For all other declarations int,FILE etc

#include "EIMTomo/common/EIMTomoTypes.h"

#include <time.h>

clock_t startm, stopm;
#define START if ( (startm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define STOP if ( (stopm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define PRINTTIME printf( "%6.3f seconds used by the processor.\n", ((double)stopm-startm)/CLOCKS_PER_SEC);

#define X_STRETCH 1
#define Z_STRETCH 2

typedef double DATA_TYPE;

#define PI 4*atan(1)

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
		DATA_TYPE delta_r;//Distance between successive measurements along x
		DATA_TYPE delta_t;//Distance between successive measurements along y
		DATA_TYPE ***counts;//The measured images should be stored in this once read from the input file. It will be a Ny X (Nz X Nx)
		DATA_TYPE *angles;//Holds the angles through which the object is tilted
		DATA_TYPE R0,RMax;
		DATA_TYPE T0,TMax;
		DATA_TYPE InitialGain,InitialOffset;//Initial scale and offset of the sinogram data 		
	};
	
	typedef struct _sinogram Sino;
	
	struct _geometry 
	{
		//User Input
		DATA_TYPE LengthZ;//This is the sample thickness
		DATA_TYPE delta_xz;//Voxel size in the x-z plane (assuming square shaped voxels in the x-z plane)
		DATA_TYPE delta_xy;//Voxel size in the x-y plane
		DATA_TYPE ***Object;//Holds the volume to be reconstructed 
		//Computed From User Input
		DATA_TYPE LengthX;//sinogram.N_x * delta_r; 
		DATA_TYPE LengthY;//sinogram.N_y * delta_t
		uint16_t N_x;//Number of voxels in x direction 
		uint16_t N_z;//Number of voxels in z direction
		uint16_t N_y;//Number of measurements in y direction
		//Coordinates of the left corner of the x-z object 
		DATA_TYPE x0;// -LengthX/2
		DATA_TYPE z0;// -LengthZ/2
		DATA_TYPE y0;//-LengthY/2
		uint16_t StartSlice,EndSlice;//Indicates Region of the object to reconstruct		
	};
	
	typedef struct _geometry Geom;
	
	struct _command_line_inputs
	{
		char* ParamFile;
		char* SinoFile;
		char* InitialRecon;
		char* OutputFile;
		char* InitialParameters;//a file containing initial gains and offsets 
		//This is read from the paramter file
		int16_t NumIter;
		DATA_TYPE SigmaX;
		DATA_TYPE p;
	};
	typedef struct _command_line_inputs CommandLineInputs;
	
	int CI_ParseInput(int ,char**,CommandLineInputs*);
	void CI_ReadParameterFile(FILE* ,CommandLineInputs* ,Sino* ,Geom*);
	void CI_InitializeSinoParameters(Sino *,CommandLineInputs*);
	void CI_InitializeGeomParameters(Sino* ,Geom* ,CommandLineInputs*);
	
	
#ifdef __cplusplus
}
#endif

#endif //ComputationInputs
