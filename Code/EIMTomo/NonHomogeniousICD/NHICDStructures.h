/*
 * NHICDStructures.h
 *
 *  Created on: Nov 11, 2011
 *      Author: mjackson
 */

#ifndef NHICDSTRUCTURES_H_
#define NHICDSTRUCTURES_H_

#include "EIMTomo/EIMTomo.h"

typedef double DATA_TYPE;


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

  typedef struct _sinogram Sino;

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

  typedef struct _geometry Geom;

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
  typedef struct _command_line_inputs CommandLineInputs;


  //Structure to store a single column(A_i) of the A-matrix
  typedef struct
  {
    double* values;//Store the non zero entries
    uint32_t count;//The number of non zero values present in the column
    uint32_t *index;//This maps each value to its location in the column. The entries in this can vary from 0 to Sinogram.N_x Sinogram.N_theta-1
  }AMatrixCol;


#endif /* NHICDSTRUCTURES_H_ */
