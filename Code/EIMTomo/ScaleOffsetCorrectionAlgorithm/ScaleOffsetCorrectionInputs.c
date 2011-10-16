#include "ScaleOffsetCorrectionInputs.h"
#include "EIMTomo/common/allocate.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


extern int optind;
extern char *optarg;

//#define START_SLICE 0
//#define END_SLICE 348

#ifdef __cplusplus
extern "C" {
#endif

int CI_ParseInput(int argc,char **argv,CommandLineInputs* Input)
{
	int j,error;
	error=0;
	while ((j = getopt(argc, argv, ":p:s:i:o:n:l:m:")) != -1) //p is for parameter file , s for sinogram file,i for Initial Recon Data File name,o is for Output File name
	{
		switch (j)
		{
			case 'i':Input->InitialRecon = optarg;break;
			case 'o':Input->OutputFile = optarg;break;
			case 'p':Input->ParamFile = optarg;break;
			case 's':Input->SinoFile = optarg;break;
			case 'n':Input->NumIter = atoi(optarg);break; //n = number of iterations
			case 'l':Input->SigmaX = (DATA_TYPE)atof(optarg);break; //l - lambda
			case 'm':Input->p = (DATA_TYPE)atof(optarg);break;//p = Markov Radom Field Parameter
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

		if(strcmp("SinoN_r",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_r=atoi(temp);
			printf("Sino.N_r=%d\n",Sinogram->N_r);
		}
		else if(strcmp("SinoN_t",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_t=atoi(temp);
			printf("Sino.N_t=%d\n",Sinogram->N_t);
		}
		else if(strcmp("SinoN_theta",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->N_theta=atoi(temp);
			printf("Sino.N_theta=%d\n",Sinogram->N_theta);
			Sinogram->angles = (DATA_TYPE*)get_spc(Sinogram->N_theta,sizeof(DATA_TYPE));
		}
		else if(strcmp("SinoDelta_r",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->delta_r=(DATA_TYPE)atof(temp);
			printf("Sino.delta_r=%lf\n",Sinogram->delta_r);
			Sinogram->R0 = -(Sinogram->N_r*Sinogram->delta_r)/2;
			Sinogram->RMax = (Sinogram->N_r*Sinogram->delta_r)/2;
		}
		else if (strcmp("SinoDelta_t",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Sinogram->delta_t=(DATA_TYPE)atof(temp);
			printf("Sino.delta_t=%lf\n",Sinogram->delta_t);
			Sinogram->T0 =  -(Sinogram->N_t*Sinogram->delta_t)/2;
			Sinogram->TMax = (Sinogram->N_t*Sinogram->delta_t)/2;

		}
		else if (strcmp("SinoAngles",temp) == 0)
		{
			for(i=0;i<Sinogram->N_theta;i++)
			{
				fscanf(Fp,"%s",temp);
				Sinogram->angles[i]=-(DATA_TYPE)atof(temp); //The theta required here is -tilt angle
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
			Geometry->delta_xz=(DATA_TYPE)atof(temp);
			printf("Geom.delta_xz=%lf\n",Geometry->delta_xz);
		}

		else if (strcmp("GeomDeltaXY",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->delta_xy=(DATA_TYPE)atof(temp);
			printf("Geom.delta_xy=%lf\n",Geometry->delta_xy);
		}
		else if (strcmp("GeomLengthZ",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->LengthZ=(DATA_TYPE)atof(temp);
			printf("Geom.LengthZ=%lf\n",Geometry->LengthZ);
		}
		else if (strcmp("StartSlice",temp) == 0)
		{
			fscanf(Fp,"%s",temp);
			Geometry->StartSlice=atoi(temp);
			printf("Geom.StartSlice=%d\n",Geometry->StartSlice);
		}
		 else if (strcmp("EndSlice",temp) == 0)
		 {
		 fscanf(Fp,"%s",temp);
		 Geometry->EndSlice=atoi(temp);
		 printf("Geom.EndSlice=%d\n",Geometry->EndSlice);
		 }
		 else if (strcmp("InitialGain",temp) == 0)
		 {
		 fscanf(Fp,"%s",temp);
		 Sinogram->InitialGain=atof(temp);
		 printf("Sinogram.InitialGain=%lf\n",Sinogram->InitialGain);
		 }
		 else if (strcmp("InitialOffset",temp) == 0)
		 {
		 fscanf(Fp,"%s",temp);
		 Sinogram->InitialOffset=atof(temp);
		 printf("Sinogram.InitialOffset=%lf\n",Sinogram->InitialOffset);
		 }
		 

	}

}

void CI_InitializeSinoParameters(Sino* Sinogram,CommandLineInputs* ParsedInput)
{
  int16_t i,j,k;
  FILE* Fp;
  double *buffer=(double*)get_spc(1,sizeof(double));
  DATA_TYPE sum=0;

  //Allocate a 3-D matrix to store the singoram in the form of a N_y X N_theta X N_x  matrix
//  Sinogram->counts=(DATA_TYPE***)get_3D(Sinogram->N_t,Sinogram->N_theta,Sinogram->N_r, sizeof(DATA_TYPE));
	 Sinogram->counts=(DATA_TYPE***)get_3D(Sinogram->N_theta,Sinogram->N_r,Sinogram->N_t, sizeof(DATA_TYPE));
  //Read data into this matrix
  //TODO: clarify this ! Super important !
 Fp=fopen(ParsedInput->SinoFile,"r");
	/*
  for(i=0;i<Sinogram->N_t;i++)
    for(j=0;j<Sinogram->N_r;j++)
      for(k=0;k<Sinogram->N_theta;k++)
  {
    fread (buffer,sizeof(DATA_TYPE),1,Fp);
    Sinogram->counts[i][k][j]=*buffer;
  }
*/
	for(i=0;i<Sinogram->N_t;i++)
		for(j=0;j<Sinogram->N_r;j++)
			for(k=0;k<Sinogram->N_theta;k++)
			{
				fread (buffer,sizeof(double), 1, Fp);
				Sinogram->counts[k][j][i] = (DATA_TYPE)(*buffer);
			}
      //check sum calculation
  for(i=0;i<Sinogram->N_theta;i++)
  {
    sum=0;
    for(j=0;j<Sinogram->N_r;j++)
      for(k=0;k<Sinogram->N_t;k++)
	sum+=Sinogram->counts[i][j][k];
    printf("Sinogram Checksum %f\n",sum);
  }
  //end ofcheck sum

  fclose(Fp);

}

void CI_InitializeGeomParameters(Sino* Sinogram,Geom* Geometry,CommandLineInputs* ParsedInput)
{
  FILE* Fp;
  uint16_t i,j,k;
  double *buffer = (double*)get_spc(1,sizeof(double));
  DATA_TYPE sum=0;//check sum TODO delete this later
  Geometry->LengthZ*=Z_STRETCH; 
	Geometry->LengthX = ((Sinogram->N_r * Sinogram->delta_r)/cos(Sinogram->angles[0]*PI/180)) + Geometry->LengthZ*tan(Sinogram->angles[0]*PI/180);
  Geometry->LengthY = (Geometry->EndSlice- Geometry->StartSlice)*Geometry->delta_xy;
  
  Geometry->N_x = ceil(Geometry->LengthX/Geometry->delta_xz);//Number of voxels in x direction
  Geometry->N_z = ceil(Geometry->LengthZ/Geometry->delta_xz);//Number of voxels in z direction
  Geometry->N_y = ceil(Geometry->LengthY/Geometry->delta_xy);//Number of measurements in y direction

	printf("Geometry->Nz=%d\n",Geometry->N_z);
	printf("Geometry->Nx=%d\n",Geometry->N_x);
	printf("Geometry->Ny=%d\n",Geometry->N_y);
	
//  Geometry->Object = (DATA_TYPE ***)get_3D(Geometry->N_y,Geometry->N_z,Geometry->N_x,sizeof(DATA_TYPE));//Allocate space for the 3-D object
  Geometry->Object = (DATA_TYPE ***)get_3D(Geometry->N_z,Geometry->N_x,Geometry->N_y,sizeof(DATA_TYPE));//Allocate space for the 3-D object
//Coordinates of the left corner of the x-z object
  Geometry->x0 = -Geometry->LengthX/2;
  Geometry->z0 = -Geometry->LengthZ/2;
  Geometry->y0 = -(Sinogram->N_t * Sinogram->delta_t)/2 + Geometry->StartSlice*Geometry->delta_xy;

  //Read the Initial Reconstruction data into a 3-D matrix
  Fp=fopen(ParsedInput->InitialRecon,"r");
/*  for(i=0;i<Geometry->N_y;i++)
    for(j=0;j<Geometry->N_x;j++)
      for(k=0;k<Geometry->N_z;k++)
  {
    fread (buffer,sizeof(DATA_TYPE),1,Fp);
    Geometry->Object[i][k][j]=*buffer;
//	printf("%f\n",Geometry->Object[i][j][k]);
  }*/

	for (i = 0; i < Geometry->N_y; i++) {
		for (j = 0; j < Geometry->N_x; j++) {
			for (k = 0; k < Geometry->N_z; k++) {
				fread(buffer, sizeof(double), 1, Fp);
				Geometry->Object[k][j][i] =  0;//(DATA_TYPE)(*buffer);
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
    printf("Geometry check sum %f\n",sum);
  }
      //End of check sum

  fclose(Fp);

}


#ifdef __cplusplus
}
#endif
