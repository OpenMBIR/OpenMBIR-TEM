/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Singanallur Venkatakrishnan (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Pudue
 * Univeristy, BlueQuartz Software nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SCALEOFFSETMOTIONSTRUCTURES_H_
#define SCALEOFFSETMOTIONSTRUCTURES_H_

#include <string>
#include <vector>

typedef double Real_t;

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/TomoArray.hpp"

typedef TomoArray<uint8_t, uint8_t*, 2> UInt8Image_t;
typedef TomoArray<uint8_t, uint8_t*, 1> UInt8ArrayType;


//typedef TomoArray<int32_t, int32_t***, 3> Int32VolumeType;

typedef TomoArray<int32_t, int32_t*, 1> Int32ArrayType;

typedef TomoArray<Real_t, Real_t*, 3> RealVolumeType;
typedef TomoArray<Real_t, Real_t*, 2> RealImage_t;
typedef TomoArray<Real_t, Real_t*, 1> RealArrayType;


namespace SOC {
  enum TiltSelection {
    A_Tilt,
    B_Tilt
  };
}

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
    uint16_t N_t;//Number of measurements in y direction
    uint16_t N_theta;//Number of angles
    Real_t delta_r;//Distance between successive measurements along x
    Real_t delta_t;//Distance between successive measurements along y
    RealVolumeType::Pointer counts;//The measured images should be stored in this once read from the input file. It will be a Ny X (Nz X Nx)
    std::vector<Real_t> angles;//Holds the angles through which the object is tilted
    Real_t R0,RMax;
    Real_t T0,TMax;
    Real_t targetGain;//,InitialOffset;//Initial scale and offset of the sinogram data
    bool BF_Flag;

    RealArrayType::Pointer InitialGain;//Reads in the initial value for the gain for each view
    RealArrayType::Pointer InitialOffset;
    RealArrayType::Pointer InitialVariance;

  } Sinogram;

  typedef boost::shared_ptr<Sinogram> SinogramPtr;

  typedef struct
  {
    RealVolumeType::Pointer Object;//Holds the volume to be reconstructed
    //Computed From User Input
    Real_t LengthX;//sinogram.N_x * delta_r;
    Real_t LengthY;//sinogram.N_y * delta_t
    uint16_t N_x;//Number of voxels in x direction
    uint16_t N_z;//Number of voxels in z direction
    uint16_t N_y;//Number of measurements in y direction
    //Coordinates of the left corner of the x-z object
    Real_t x0;// -LengthX/2
    Real_t z0;// -LengthZ/2
    Real_t y0;//-LengthY/2
  } Geometry;

  typedef boost::shared_ptr<Geometry> GeometryPtr;

  typedef struct
  {
    uint16_t NumIter;
    uint16_t NumOuterIter;
    Real_t SigmaX;
    Real_t p;
    Real_t StopThreshold;
    uint8_t   InterpFlag;

    bool useSubvolume;
    uint16_t xStart;
    uint16_t xEnd;
    uint16_t yStart;
    uint16_t yEnd;
    uint16_t zStart;
    uint16_t zEnd;

    SOC::TiltSelection tiltSelection;

    uint16_t fileXSize; // Size in voxels of the complete width of the image from the file
    uint16_t fileYSize; // Size in voxels of the complete height of the image from the file
    uint16_t fileZSize; // Number of tilts in the series from the file.

    Real_t LengthZ;//This is the sample thickness
    Real_t delta_xz;//Voxel size in the x-z plane (assuming square shaped voxels in the x-z plane)
    Real_t delta_xy;//Voxel size in the x-y plane

  	bool extendObject; //In case the sinogram data corresponding to voxels outside it
    Real_t interpolateFactor;
    Real_t targetGain;
    bool useDefaultOffset;
    Real_t defaultOffset;
    Real_t defaultInitialRecon;
    Real_t defaultVariance;


    /* These are input files */
    std::string sinoFile; /* .mrc formatted files are accepted currently */
    std::string initialReconFile;

    std::string gainsInputFile;
    std::string offsetsInputFile;
    std::string varianceInputFile;

    /* These are output related files and parameters */
    std::string tempDir; // Output directory
    std::string reconstructedOutputFile;
    std::string gainsOutputFile;
    std::string offsetsOutputFile;
    std::string varianceOutputFile;


    std::vector<uint8_t> excludedViews;// Indices of views to exclude from reconstruction
    std::vector<int> goodViews; // Contains the indices of the views to use for reconstruction
  } TomoInputs;
  typedef boost::shared_ptr<TomoInputs> TomoInputsPtr;


  typedef struct
  {
      double X_SHRINK_FACTOR;
      unsigned int X_STRETCH;
      unsigned int Z_STRETCH;
      unsigned int DETECTOR_RESPONSE_BINS;
      unsigned int PROFILE_RESOLUTION;
      unsigned int BEAM_RESOLUTION;
      unsigned int AREA_WEIGHTED;
      unsigned int THRESHOLD_REDUCTION_FACTOR;
      unsigned int ZERO_SKIPPING;
  } AdvancedParameters ;
  typedef boost::shared_ptr<AdvancedParameters> AdvancedParametersPtr;

  //Structure to store a single column(A_i) of the A-matrix
  typedef struct
  {
      Real_t* values; //Store the non zero entries
      uint32_t count; //The number of non zero values present in the column
      uint32_t *index; //This maps each value to its location in the column. The entries in this can vary from 0 to Sinogram.N_x Sinogram.N_theta-1
  } AMatrixCol;

  typedef struct
  {

    RealArrayType::Pointer I_0; //Gains
    RealArrayType::Pointer mu; //Offset
    RealArrayType::Pointer alpha;//Noise variance refinement factor
  } ScaleOffsetParams;
  typedef boost::shared_ptr<ScaleOffsetParams> ScaleOffsetParamsPtr;


#endif /* SCALEOFFSETMOTIONSTRUCTURES_H_ */
