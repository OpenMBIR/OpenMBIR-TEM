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

#ifndef _ReconstructionStructures_H_
#define _ReconstructionStructures_H_

#include <string>
#include <vector>


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/TomoArray.hpp"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"


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
  Real_t R0;
  Real_t RMax;
  Real_t T0;
  Real_t TMax;
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
  int16_t  NumThreads;
  Real_t SigmaX;
  Real_t p;
  Real_t StopThreshold;
  uint8_t   InterpFlag;
  bool extendObject; //In case the sinogram data corresponding to voxels outside it
  bool verbose;
  bool veryVerbose;
  bool useSubvolume;
  uint16_t xStart;
  uint16_t xEnd;
  uint16_t yStart;
  uint16_t yEnd;
  uint16_t zStart;
  uint16_t zEnd;

  uint16_t fileXSize; // Size in voxels of the complete width of the image from the file
  uint16_t fileYSize; // Size in voxels of the complete height of the image from the file
  uint16_t fileZSize; // Number of tilts in the series from the file.

  Real_t LengthZ;//This is the sample thickness
  Real_t delta_xz;//Voxel size in the x-z plane (assuming square shaped voxels in the x-z plane)
  Real_t delta_xy;//Voxel size in the x-y plane

  Real_t interpolateFactor;
  Real_t defaultInitialRecon;

  std::vector<float> tilts;

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
  std::string vtkOutputFile;
  std::string mrcOutputFile;
  std::string avizoOutputFile;
  std::string braggSelectorFile;
  std::vector<std::string> tempFiles;

  std::vector<uint8_t> excludedViews;// Indices of views to exclude from reconstruction
  std::vector<int> goodViews; // Contains the indices of the views to use for reconstruction

  /* These are Mainly for the HAADF reconstructions */
  Real_t targetGain;
  bool useDefaultOffset;
  Real_t defaultOffset;
  Real_t defaultVariance;

} TomoInputs;
typedef boost::shared_ptr<TomoInputs> TomoInputsPtr;


typedef struct
{
  double X_SHRINK_FACTOR; /* This is defaulted to 0.6 */
  unsigned int X_STRETCH; /* This is defaulted to 1 */
  unsigned int Z_STRETCH; /* This is defaulted to 2 */
  uint16_t DETECTOR_RESPONSE_BINS; /* This is defaulted to 64 */
  int32_t PROFILE_RESOLUTION; /* This is defaulted to 1536 */
  unsigned int BEAM_RESOLUTION; /* This is defaulted to 512 */
  unsigned int AREA_WEIGHTED; /* This is defaulted to 1 */
  unsigned int THRESHOLD_REDUCTION_FACTOR; /* This is defaulted to 1 */
  unsigned int JOINT_ESTIMATION; /* It need not be on for all cases. In case
                                    the Bright Field image file is present (very unlikely) + the offset
                                    is known the user might want to Turn this OFF. All
                                    other cases it needs be switched ON. */
  unsigned int ZERO_SKIPPING; /* This will always be ON in the end product I think. */
  unsigned int NOISE_ESTIMATION; /* This is a parameter that the user MAY or MAY NOT
                                  want turned ON. It is ON by default */
  unsigned int NOISE_MODEL; /* This is a parameter that the user MAY or MAY NOT
                want turned ON. It is ON by default */
  unsigned int ESTIMATE_PRIOR;
} AdvancedParameters;
typedef boost::shared_ptr<AdvancedParameters> AdvancedParametersPtr;

#endif /* _ReconstructionStructures_H_ */
