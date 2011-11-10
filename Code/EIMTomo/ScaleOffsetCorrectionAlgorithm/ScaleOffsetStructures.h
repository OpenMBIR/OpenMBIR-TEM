

#ifndef SCALEOFFSETSTRUCTURES_H_
#define SCALEOFFSETSTRUCTURES_H_




typedef double DATA_TYPE;

#include "EIMTomo/EIMTomo.h"


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

  typedef struct
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
    DATA_TYPE TargetGain;//,InitialOffset;//Initial scale and offset of the sinogram data
    uint8_t *ViewMask;//Which views to keep and which to reject
    uint16_t N_tStart,N_tEnd,N_rStart,N_rEnd;//Which region of the sinogram to keep in the r and t directions
    DATA_TYPE *InitialGain;//Reads in the initial value for the gain for each view
    DATA_TYPE *InitialOffset;
  } Sino;



  typedef struct
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
  } Geom;


  typedef struct
  {
    char* ParamFile;
    char* SinoFile;
    char* InitialRecon;
    char* OutputFile;
    char* InitialParameters;//a file containing initial gains and offsets
    //This is read from the paramter file
    int16_t NumIter;
    uint16_t NumOuterIter;
    DATA_TYPE SigmaX;
    DATA_TYPE p;
  } CommandLineInputs;

  //Structure to store a single column(A_i) of the A-matrix
  typedef struct
  {
      DATA_TYPE* values; //Store the non zero entries
      uint32_t count; //The number of non zero values present in the column
      uint32_t *index; //This maps each value to its location in the column. The entries in this can vary from 0 to Sinogram.N_x Sinogram.N_theta-1
  } AMatrixCol;

  typedef struct
  {
      DATA_TYPE *I_0; //Scale
      DATA_TYPE *mu; //Offset
  } ScaleOffsetParams;



#endif /* SCALEOFFSETSTRUCTURES_H_ */
