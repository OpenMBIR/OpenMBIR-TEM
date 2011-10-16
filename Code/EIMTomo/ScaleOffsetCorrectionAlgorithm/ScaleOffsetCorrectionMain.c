//#include "ComputationEngine.h"
#include "ScaleOffsetCorrectionInputs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>




/*double
 CE => Computation Engine
 CI => Computation Inputs
 N_ => Number of ..
 
 */



int main(int argc,char** argv)
{
	int16_t error,i,j,k;
	FILE* Fp;
	CommandLineInputs ParsedInput;
	Sino Sinogram;
	Geom Geometry;
	DATA_TYPE *buffer=(DATA_TYPE*)get_spc(1,sizeof(DATA_TYPE));
	
	
	START;
	
	error=-1;
	error=CI_ParseInput(argc,argv,&ParsedInput);
	if(error != -1)
	{
		printf("%s %s %s %s\n",ParsedInput.ParamFile,ParsedInput.SinoFile,ParsedInput.InitialRecon,ParsedInput.OutputFile);
		//Read the paramters into the structures
		Fp = fopen(ParsedInput.ParamFile,"r");
		CI_ReadParameterFile(Fp,&ParsedInput,&Sinogram,&Geometry);
		fclose(Fp);
		//Based on the inputs , calculate the "other" variables in the structure definition
		CI_InitializeSinoParameters(&Sinogram,&ParsedInput);
		CI_InitializeGeomParameters(&Sinogram,&Geometry,&ParsedInput);
	}
			
	error=CE_MAPICDReconstruct(&Sinogram,&Geometry,&ParsedInput);
	Fp=fopen(ParsedInput.OutputFile,"w");
	
	printf("Main\n");
	printf("Final Dimensions of Object Nz=%d Nx=%d Ny=%d\n",Geometry.N_z,Geometry.N_x,Geometry.N_y);
	
	for(i = 0;i < Geometry.N_y; i++)
	{
		for(j = 0;j < Geometry.N_x; j++)
			for(k = 0;k < Geometry.N_z; k++)
			{
				buffer = &Geometry.Object[k][j][i];
				fwrite(buffer,sizeof(DATA_TYPE),1,Fp);
			}
		printf("%d\n",i);
	}
	
	fclose(Fp);
	STOP;
	PRINTTIME; 
	// if(error < 0)
	// {
	// //TODO:clean up memory here 
	// return EXIT_FAILURE;
	// }
	// //TODO:free memory
	// 
	// return EXIT_SUCCESS;
	return 0;
}

