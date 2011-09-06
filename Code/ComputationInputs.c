#include "ComputationInputs.h"
#include "allocate.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//#ifdef __linux__
#include <getopt.h>

extern int optind;
extern char *optarg;

#ifdef __cplusplus
extern "C" {
#endif

int CI_ParseInput(int argc,char **argv,CommandLineInputs* Input)
{
	int i,error;
	error=0;
	while ((i = getopt(argc, argv, ":p:s:i:o:n:l:m:")) != -1) //p is for parameter file , s for sinogram file,i for Initial Recon Data File name,o is for Output File name  
	{
		switch (i) 
		{
			case 'i':Input->InitialRecon = optarg;break;
			case 'o':Input->OutputFile = optarg;break;
			case 'p':Input->ParamFile = optarg;break;
			case 's':Input->SinoFile = optarg;break;
			case 'n':Input->NumIter = atoi(optarg);break; //n = number of iterations
			case 'l':Input->SigmaX = (double)atof(optarg);break; //l - lambda
			case 'm':Input->p = (double)atof(optarg);break;//p = Markov Radom Field Parameter
			case '?':error=-1;break;
		}
	}
	
	return error;
}

void CI_ReadParameterFile(FILE *Fp,CommandLineInputs* ParsedInput,Sino* Sinogram,Geom* Geometry)
{
	char temp[20];
	int16_t i;
	while(!feof(Fp))
	{
		fscanf(Fp,"%s",temp);
		
		if(strcmp("SinoN_x",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_x=atoi(temp);
			printf("Sino.N_x=%d\n",Sinogram->N_x);
		}
		else if(strcmp("SinoN_y",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_y=atoi(temp);
			printf("Sino.N_y=%d\n",Sinogram->N_y);
		}
		else if(strcmp("SinoN_theta",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_theta=atoi(temp);
			printf("Sino.N_theta=%d\n",Sinogram->N_theta);
			Sinogram->angles = (double*)get_spc(Sinogram->N_theta,sizeof(double));
		}
		else if(strcmp("SinoDelta_x",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->delta_x=(double)atof(temp);
			printf("Sino.delta_x=%lf\n",Sinogram->delta_x);
			Sinogram->R0 = -(Sinogram->N_x*Sinogram->delta_x)/2;
			Sinogram->RMax = (Sinogram->N_x*Sinogram->delta_x)/2;
		}
		else if (strcmp("SinoDelta_y",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->delta_y=(double)atof(temp);
			printf("Sino.delta_y=%lf\n",Sinogram->delta_y);
			Sinogram->T0 =  -(Sinogram->N_y*Sinogram->delta_y)/2;
			Sinogram->TMax = (Sinogram->N_y*Sinogram->delta_y)/2; 
			
		}
		else if (strcmp("SinoAngles",temp) == 0)
		{
			for(i=0;i<Sinogram->N_theta;i++)
			{
				fscanf(Fp,"%s",temp);
				Sinogram->angles[i]=-(double)atof(temp); //The theta required here is -tilt angle
				printf("%lf\n",Sinogram->angles[i]);
			}
		}
	/*	else if (strcmp("Iter",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			ParsedInput->NumIter=atoi(temp);
			printf("Params.NumIter=%d\n",ParsedInput->NumIter);
		}*/
		
		else if (strcmp("GeomDeltaXZ",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->delta_xz=(double)atof(temp);
			printf("Geom.delta_xz=%lf\n",Geometry->delta_xz);
		}    
		
		else if (strcmp("GeomDeltaXY",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->delta_xy=(double)atof(temp);
			printf("Geom.delta_xy=%lf\n",Geometry->delta_xy);
		}
		else if (strcmp("GeomLengthZ",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->LengthZ=(double)atof(temp);
			printf("Geom.LengthZ=%lf\n",Geometry->LengthZ);
		}        
		/*else if (strcmp("Scaling",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			ParsedInput->scaling=atof(temp);
			printf("Params.scaling=%lf\n",ParsedInput->scaling);
		}*/
		
	}
	
}

void CI_InitializeSinoParameters(Sino* Sinogram,CommandLineInputs* ParsedInput)
{
  int16_t i,j,k;
  FILE* Fp;
  double *buffer=(double*)get_spc(1,sizeof(float));
  double sum=0;
  
  //Allocate a 3-D matrix to store the singoram in the form of a N_y X N_theta X N_x  matrix 
//  Sinogram->counts=(double***)get_3D(Sinogram->N_y,Sinogram->N_theta,Sinogram->N_x, sizeof(double));
	 Sinogram->counts=(double***)get_3D(Sinogram->N_theta,Sinogram->N_x,Sinogram->N_y, sizeof(double));
  //Read data into this matrix
  //TODO: clarify this ! Super important ! 
 Fp=fopen(ParsedInput->SinoFile,"r");
	/*
  for(i=0;i<Sinogram->N_y;i++)
    for(j=0;j<Sinogram->N_x;j++)
      for(k=0;k<Sinogram->N_theta;k++)
  {	 
    fread (buffer,sizeof(double),1,Fp);
    Sinogram->counts[i][k][j]=*buffer;
  }
*/
	for(i=0;i<Sinogram->N_y;i++)
		for(j=0;j<Sinogram->N_x;j++)
			for(k=0;k<Sinogram->N_theta;k++)
			{	 
				fread (buffer,sizeof(double), 1, Fp);
				Sinogram->counts[k][j][i] = *buffer;
			}
      //check sum calculation
  for(i=0;i<Sinogram->N_theta;i++)
  {
    sum=0;
    for(j=0;j<Sinogram->N_x;j++)
      for(k=0;k<Sinogram->N_y;k++)
	sum+=Sinogram->counts[i][j][k];
    printf("Sinogram Checksum %lf\n",sum);
  }
  //end ofcheck sum
	    
  fclose(Fp);
  
}

void CI_InitializeGeomParameters(Sino* Sinogram,Geom* Geometry,CommandLineInputs* ParsedInput)
{
  FILE* Fp;
  uint16_t i,j,k;
  double *buffer = (double*)get_spc(1,sizeof(double));
  double sum=0;//check sum TODO delete this later
  Geometry->LengthX = Sinogram->N_x * Sinogram->delta_x;//sinogram.N_x * delta_x; 
  Geometry->LengthY = Sinogram->N_y * Sinogram->delta_y;//sinogram.N_y * delta_y
  Geometry->N_x = ceil(Geometry->LengthX/Geometry->delta_xz);//Number of voxels in x direction 
  Geometry->N_z = ceil(Geometry->LengthZ/Geometry->delta_xz);//Number of voxels in z direction
  Geometry->N_y = ceil(Geometry->LengthY/Geometry->delta_xy);//Number of measurements in y direction
//  Geometry->Object = (double ***)get_3D(Geometry->N_y,Geometry->N_z,Geometry->N_x,sizeof(double));//Allocate space for the 3-D object
  Geometry->Object = (double ***)get_3D(Geometry->N_z,Geometry->N_x,Geometry->N_y,sizeof(double));//Allocate space for the 3-D object
//Coordinates of the left corner of the x-z object 
  Geometry->x0 = -Geometry->LengthX/2;
  Geometry->z0 = -Geometry->LengthZ/2;
	Geometry->y0 = -Geometry->LengthY/2;
  
  //Read the Initial Reconstruction data into a 3-D matrix
  Fp=fopen(ParsedInput->InitialRecon,"r");
/*  for(i=0;i<Geometry->N_y;i++)
    for(j=0;j<Geometry->N_x;j++)
      for(k=0;k<Geometry->N_z;k++)
  {	 
    fread (buffer,sizeof(double),1,Fp);
    Geometry->Object[i][k][j]=*buffer;
//	printf("%f\n",Geometry->Object[i][j][k]);
  }*/
	
	for (i = 0; i < Geometry->N_y; i++) {
		for (j = 0; j < Geometry->N_x; j++) {
			for (k = 0; k < Geometry->N_z; k++) {
				fread(buffer, sizeof(double), 1, Fp);
				Geometry->Object[k][j][i] = *buffer;
			}
		}
	}

      //Doing a check sum to verify with matlab
  
  for(i=0;i<Geometry->N_y;i++)
  {
    sum=0;
    for(j=0;j<Geometry->N_x;j++)
      for(k=0;k<Geometry->N_z;k++)
    {
      sum+=Geometry->Object[k][j][i];
    }
    printf("Geometry check sum %lf\n",sum);
  }
      //End of check sum
     
  fclose(Fp);
 
}


#ifdef __cplusplus
}
#endif
