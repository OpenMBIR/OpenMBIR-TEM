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
	while ((j = getopt(argc, argv, ":p:s:i:o:n:l:m:g:O:")) != -1) //p is for parameter file , s for sinogram file,i for Initial Recon Data File name,o is for Output File name
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
			case 'g':Input->InitialParameters = optarg;break;//the values of the initial gain and offset
			case 'O':Input->NumOuterIter = atoi(optarg);break;//Number of Outer Iterations for the joint estimation
			case '?':error=-1;break;
		}
	}

	return error;
}

void CI_ReadParameterFile(FILE *Fp,CommandLineInputs* ParsedInput,Sino* Sinogram,Geom* Geometry)
{
	char temp[20];
	int16_t i;
	uint16_t view_count;
	DATA_TYPE tempvariable=0;
	DATA_TYPE *MaskedViewAngles;
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
			Sinogram->ViewMask = (DATA_TYPE*)get_spc(Sinogram->N_theta, sizeof(uint8_t));

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
			view_count=0;
			for(i=0;i<Sinogram->N_theta;i++)
			{
				fscanf(Fp,"%s",temp);
				Sinogram->angles[i]=-(DATA_TYPE)atof(temp); //The theta required here is -tilt angle
				printf("%lf\n",Sinogram->angles[i]);
				fscanf(Fp,"%s",temp);
				Sinogram->ViewMask[i]=atoi(temp); //A mask to select or reject the view i
				if(Sinogram->ViewMask[i] == 1)
					view_count++;

			//	printf("Mask %d\n",Sinogram->ViewMask[i]);
			}
			printf("Number Of Views %d\n",view_count);
			MaskedViewAngles = (DATA_TYPE*)get_spc(view_count, sizeof(DATA_TYPE));

			view_count=0;
			for(i=0;i<Sinogram->N_theta;i++)
			{
				if (Sinogram->ViewMask[i] ==1)
				{
					MaskedViewAngles[view_count++]=Sinogram->angles[i];
				}
			}

			free(Sinogram->angles);//clear the initial set of view angles to replace it with the masked set

		   Sinogram->angles = (DATA_TYPE*)get_spc(view_count,sizeof(DATA_TYPE));

			for (i=0; i < view_count;i++)
			{
				Sinogram->angles[i] = MaskedViewAngles[i];
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

		/*
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
		 */

		 else if (strcmp("TargetGain",temp) == 0)
		 {
		 fscanf(Fp,"%s",temp);
		 Sinogram->TargetGain=atof(temp);
		 printf("Sinogram.TargetGain=%lf\n",Sinogram->TargetGain);
		 }
		/*
		 else if (strcmp("InitialOffset",temp) == 0)
		 {
		 fscanf(Fp,"%s",temp);
		 Sinogram->InitialOffset=atof(temp);
		 printf("Sinogram.InitialOffset=%lf\n",Sinogram->InitialOffset);
		 }*/
		 else if (strcmp("N_tStart",temp) == 0)
		 {
			 fscanf(Fp,"%s",temp);
			 Sinogram->N_tStart=atoi(temp);
			 printf("Sino.N_tStartSlice=%d\n",Sinogram->N_tStart);
		 }
		 else if (strcmp("N_tEnd",temp) == 0)
		 {
			 fscanf(Fp,"%s",temp);
			 Sinogram->N_tEnd=atoi(temp);
			 printf("Sino.N_tEndSlice=%d\n",Sinogram->N_tEnd);
		 }
		 else if (strcmp("N_rStart",temp) == 0)
		 {
			 fscanf(Fp,"%s",temp);
			 Sinogram->N_rStart=atoi(temp);
			 printf("Sino.N_rStartSlice=%d\n",Sinogram->N_rStart);
		 }
		 else if (strcmp("N_rEnd",temp) == 0)
		 {
			 fscanf(Fp,"%s",temp);
			 Sinogram->N_rEnd=atoi(temp);
			 printf("Sino.N_rEndSlice=%d\n",Sinogram->N_rEnd);
		 }

	}

}

void CI_InitializeSinoParameters(Sino* Sinogram,CommandLineInputs* ParsedInput)
{
  int16_t i,j,k;
	uint16_t view_count=0,TotalNumMaskedViews;
  FILE* Fp;
  double *buffer=(double*)get_spc(1,sizeof(double));
  DATA_TYPE sum=0;

	for(i=0; i < Sinogram->N_theta;i++)
		if(Sinogram->ViewMask[i] == 1)
			view_count++;
	TotalNumMaskedViews=view_count;

  //Allocate a 3-D matrix to store the singoram in the form of a N_y X N_theta X N_x  matrix
//  Sinogram->counts=(DATA_TYPE***)get_3D(Sinogram->N_t,Sinogram->N_theta,Sinogram->N_r, sizeof(DATA_TYPE));
	 Sinogram->counts=(DATA_TYPE***)get_3D(TotalNumMaskedViews,Sinogram->N_rEnd-Sinogram->N_rStart+1,Sinogram->N_tEnd-Sinogram->N_rStart+1, sizeof(DATA_TYPE));
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
		{
			view_count=0;
			for(k=0;k<Sinogram->N_theta;k++)
			{
				fread (buffer,sizeof(double), 1, Fp);
				if(Sinogram->ViewMask[k] == 1 && j >= Sinogram->N_rStart && j <= Sinogram->N_rEnd && i >= Sinogram->N_tStart && i <= Sinogram->N_tEnd)
				Sinogram->counts[view_count++][j-Sinogram->N_rStart][i-Sinogram->N_tStart] = (DATA_TYPE)(*buffer);
			}


		}
	fclose(Fp);

	//The normalization and offset parameters for the views
	Sinogram->InitialGain=(DATA_TYPE*)get_spc(TotalNumMaskedViews, sizeof(DATA_TYPE));
	Sinogram->InitialOffset=(DATA_TYPE*)get_spc(TotalNumMaskedViews, sizeof(DATA_TYPE));

	Fp=fopen(ParsedInput->InitialParameters,"r");//This file contains the Initial unscatterd counts and background scatter for each view

	view_count=0;
	for( i = 0; i < Sinogram->N_theta; i++)
	{
		fread(buffer, sizeof(double), 1, Fp);
		if(Sinogram->ViewMask[i] == 1)
			Sinogram->InitialGain[view_count++]=(DATA_TYPE)(*buffer);
	}
	view_count=0;
	for( i = 0; i < Sinogram->N_theta; i++)
	{
		fread(buffer, sizeof(double), 1, Fp);
		if(Sinogram->ViewMask[i] == 1)
			Sinogram->InitialOffset[view_count++]=(DATA_TYPE)(*buffer);
	}




	Sinogram->N_theta = TotalNumMaskedViews;
	Sinogram->N_r = (Sinogram->N_rEnd-Sinogram->N_rStart+1);
	Sinogram->N_t = (Sinogram->N_tEnd-Sinogram->N_tStart+1);
	Sinogram->R0 = -(Sinogram->N_r*Sinogram->delta_r)/2;
	Sinogram->RMax = (Sinogram->N_r*Sinogram->delta_r)/2;
	Sinogram->T0 =  -(Sinogram->N_t*Sinogram->delta_t)/2;
	Sinogram->TMax = (Sinogram->N_t*Sinogram->delta_t)/2;


	printf("Size of the Masked Sinogram N_r =%d N_t = %d N_theta=%d\n",Sinogram->N_r,Sinogram->N_t,Sinogram->N_theta);

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

	/*

//This function takes the original sinogram and creates a masked version of it based on the ROI parameters . In the final version of the code all this should be incorporated into
//FILE read instead of pruning it in the program after reading the whole input file
void CI_MaskSinogram(Sino* OriginalSinogram,Sino* NewSinogram)
{
	uint16_t i_theta,i_r,i_t,view_count=0;
	NewSinogram->N_r = (OriginalSinogram->N_rEnd-OriginalSinogram->N_rStart+1);
	NewSinogram->N_t = (OriginalSinogram->N_tEnd-OriginalSinogram->N_tStart+1);

	for(i_theta = 0;i_theta < OriginalSinogram->N_theta;i_theta++)
		view_count+=OriginalSinogram->ViewMask[i_theta];

	NewSinogram->N_theta = view_count;
	NewSinogram->angles=(DATA_TYPE*)get_spc(NewSinogram->N_theta, sizeof(DATA_TYPE));
	NewSinogram->counts=(DATA_TYPE***)get_3D(NewSinogram->N_theta,NewSinogram->N_r,NewSinogram->N_t, sizeof(DATA_TYPE));

	view_count=0;
	for(i_theta =0; i_theta < OriginalSinogram->N_theta;i_theta++)
	{
		if(OriginalSinogram->ViewMask[i_theta] == 1)//View is to be taken into account
		{
			NewSinogram->angles[view_count]=OriginalSinogram->angles[i_theta];
			for (i_r = OriginalSinogram->N_rStart; i_r <= OriginalSinogram->N_rEnd; i_r++)
				for(i_t = OriginalSinogram->N_tStart; i_t <= OriginalSinogram->N_tEnd;i_t++)
					NewSinogram[view_count][i_r-OriginalSinogram->N_rStart][i_t - OriginalSinogram->N_tStart] = OriginalSinogram->counts[i_theta][i_r][i_t];
			view_count++;
		}
	}
	NewSinogram->delta_r=OriginalSinogram->delta_r;
	NewSinogram->delta_t =OriginalSinogram->delta_t;

	NewSinogram->R0 = -(NewSinogram->N_r*NewSinogram->delta_r)/2;
	NewSinogram->RMax = (NewSinogram->N_r*NewSinogram->delta_r)/2;
	NewSinogram->T0 =  -(NewSinogram->N_t*NewSinogram->delta_t)/2;
	NewSinogram->TMax = (NewSinogram->N_t*NewSinogram->delta_t)/2;
	NewSinogram->TargetGain=OriginalSinogram->TargetGain;
	NewSinogram->InitialOffset=	OriginalSinogram->InitialOffset;

	free(OriginalSinogram->counts);//this is a huge array; eliminate it
	free(OriginalSinogram->angles);


}
	 */

void CI_InitializeGeomParameters(Sino* Sinogram,Geom* Geometry,CommandLineInputs* ParsedInput)
{
  FILE* Fp;
  uint16_t i,j,k;
  double *buffer = (double*)get_spc(1,sizeof(double));
  DATA_TYPE sum=0,max;

	//Find the maximum absolute tilt angle
	max= AbsMaxArray(Sinogram->angles, Sinogram->N_theta);

#ifndef FORWARD_PROJECT_MODE
    Geometry->LengthZ*=Z_STRETCH;

#ifdef EXTEND_OBJECT
	Geometry->LengthX = ((Sinogram->N_r * Sinogram->delta_r)/cos(max*M_PI/180)) + Geometry->LengthZ*tan(max*M_PI/180) ;
#else
	Geometry->LengthX = ((Sinogram->N_r * Sinogram->delta_r));
#endif //Extend object endif

#else
	Geometry->LengthX = ((Sinogram->N_r * Sinogram->delta_r));
#endif//Forward projector mode end if

//  Geometry->LengthY = (Geometry->EndSlice- Geometry->StartSlice)*Geometry->delta_xy;
	Geometry->LengthY = (Sinogram->N_tEnd-Sinogram->N_tStart + 1)*Sinogram->delta_t;

  Geometry->N_x = ceil(Geometry->LengthX/Geometry->delta_xz);//Number of voxels in x direction
  Geometry->N_z = ceil(Geometry->LengthZ/Geometry->delta_xz);//Number of voxels in z direction
  Geometry->N_y = floor(Geometry->LengthY/Geometry->delta_xy);//Number of measurements in y direction

	printf("Geometry->Nz=%d\n",Geometry->N_z);
	printf("Geometry->Nx=%d\n",Geometry->N_x);
    printf("Geometry->Ny=%d\n",Geometry->N_y);

//  Geometry->Object = (DATA_TYPE ***)get_3D(Geometry->N_y,Geometry->N_z,Geometry->N_x,sizeof(DATA_TYPE));//Allocate space for the 3-D object
  Geometry->Object = (DATA_TYPE ***)get_3D(Geometry->N_z,Geometry->N_x,Geometry->N_y,sizeof(DATA_TYPE));//Allocate space for the 3-D object
//Coordinates of the left corner of the x-z object
  Geometry->x0 = -Geometry->LengthX/2;
  Geometry->z0 = -Geometry->LengthZ/2;
 // Geometry->y0 = -(Sinogram->N_t * Sinogram->delta_t)/2 + Geometry->StartSlice*Geometry->delta_xy;
	Geometry->y0 = -(Geometry->LengthY)/2 ;

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

	for (i = 0; i < Geometry->N_y; i++)
	{
		for (j = 0; j < Geometry->N_x; j++)
		{
			for (k = 0; k < Geometry->N_z; k++)
			{
				if(Fp == NULL)//If no input file has been specified or if the file does not exist just set the default values to be zero
				{
				Geometry->Object[k][j][i] = 0;
				}
				else//If the iput file exists read the values
				{
				fread(buffer, sizeof(double), 1, Fp);
				Geometry->Object[k][j][i] = (DATA_TYPE)(*buffer);
				}
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

//Finds the maximum of absolute value elements in an array
DATA_TYPE AbsMaxArray(DATA_TYPE* Array ,uint16_t NumElts)
{
	uint16_t i;
	DATA_TYPE max;
	max = fabs(Array[0]);
	for(i =1; i < NumElts;i++)
		if(fabs(Array[i]) > max)
			max=fabs(Array[i]);
	return max;

}

#ifdef __cplusplus
}
#endif
