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

#include "InitialReconstructionInitializer.h"

#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Common/allocate.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionInitializer::InitialReconstructionInitializer()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InitialReconstructionInitializer::~InitialReconstructionInitializer()
{
}

// -----------------------------------------------------------------------------
//Finds the maximum of absolute value elements in an array
// -----------------------------------------------------------------------------
Real_t InitialReconstructionInitializer::absMaxArray(std::vector<Real_t> &Array)
{
  uint16_t i;
  Real_t max;
  max = fabs(Array[0]);
  for(i =1; i < Array.size();i++)
    if(fabs(Array[i]) > max)
      max=fabs(Array[i]);
  return max;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InitialReconstructionInitializer::initializeData()
{
  GeometryPtr geometry = getGeometry();

  TomoInputsPtr input = getTomoInputs();
  for (uint16_t z = 0; z < geometry->N_z; z++)
  {
    for (uint16_t x = 0; x < geometry->N_x; x++)
    {
      for (uint16_t y = 0; y < geometry->N_y; y++)
      {

        //geometry->Object->d[k][j][i] = 0;
        geometry->Object->setValue(input->defaultInitialRecon, z, x, y);
      }
    }
  }
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InitialReconstructionInitializer::execute()
{

  SinogramPtr sinogram = getSinogram();
  TomoInputsPtr input = getTomoInputs();
  GeometryPtr geometry = getGeometry();
  AdvancedParametersPtr advParams = getAdvParams();

  Real_t sum = 0, max;

#ifndef FORWARD_PROJECT_MODE
  input->delta_xz = sinogram->delta_r * input->delta_xz;
  input->delta_xy = input->delta_xz;
  //Find the maximum absolute tilt angle
  max = absMaxArray(sinogram->angles);

  //The max is going to be used to reconstruct a large area
 //However if its close to 90 this would result in a very large value - so truncate
 if(max > MBIR::Constants::k_MaxAngleStretch)
   max = MBIR::Constants::k_MaxAngleStretch;

  // Convert Max to radians
  max = max * M_PI / 180.0;
  input->LengthZ *= advParams->Z_STRETCH;
  input->LengthZ /= (input->interpolateFactor * sinogram->delta_r);
  //interpolation_factor;
  input->LengthZ = floor(input->LengthZ + 0.5) * input->interpolateFactor * sinogram->delta_r; //interpolation_factor;

  geometry->LengthX = ((sinogram->N_r * sinogram->delta_r));
  geometry->N_x = floor(geometry->LengthX / input->delta_xz); //Number of voxels in x direction

  if(1 == input->extendObject)
  {
    std::cout << "KNOWN BUG FIX NEEDED HERE IF MAX = 90 degrees" << std::endl;
    geometry->LengthX = advParams->X_SHRINK_FACTOR * ((sinogram->N_r * sinogram->delta_r) / cos(max)) + input->LengthZ * tan(max);
    geometry->LengthX /= (input->interpolateFactor * sinogram->delta_r);
    geometry->LengthX = floor(geometry->LengthX + 0.5) * input->interpolateFactor * sinogram->delta_r;
  }


#else
  geometry->LengthX = ((sinogram->N_r * sinogram->delta_r));
#endif//Forward projector mode end if
//  Geometry->LengthY = (Geometry->EndSlice- Geometry->StartSlice)*Geometry->delta_xy;
  geometry->LengthY = (input->yEnd - input->yStart + 1) * sinogram->delta_t;

  geometry->N_x = floor(geometry->LengthX / input->delta_xz); //Number of voxels in x direction
  geometry->N_z = floor(input->LengthZ / input->delta_xz); //Number of voxels in z direction
  geometry->N_y = floor(geometry->LengthY / input->delta_xy); //Number of measurements in y direction

  std::stringstream ss;

  ss << "Geometry->LengthX=" << geometry->LengthX << " nm"  << std::endl;
  ss << "Geometry->LengthY=" << geometry->LengthY << " nm"  << std::endl;
  ss << "Geometry->LengthZ=" << input->LengthZ << " nm"  << std::endl;

  ss << "Geometry->Nz=" <<  geometry->N_z << std::endl;
  ss << "Geometry->Nx=" <<  geometry->N_x << std::endl;
  ss << "Geometry->Ny=" <<  geometry->N_y << std::endl;

  size_t dims[3] =
  { geometry->N_z, geometry->N_x, geometry->N_y };
  geometry->Object = RealVolumeType::New(dims, "Geometry.Object");

  // geometry->Object = (DATA_TYPE ***)get_3D(geometry->N_z, geometry->N_x, geometry->N_y, sizeof(DATA_TYPE));//Allocate space for the 3-D object
//Coordinates of the left corner of the x-z object
  geometry->x0 = -geometry->LengthX / 2;
  geometry->z0 = -input->LengthZ / 2;
  // Geometry->y0 = -(sinogram->N_t * sinogram->delta_t)/2 + Geometry->StartSlice*Geometry->delta_xy;
  geometry->y0 = -(geometry->LengthY) / 2;

  ss << "Geometry->X0=" << geometry->x0 << std::endl;
  ss << "Geometry->Y0=" << geometry->y0 << std::endl;
  ss << "Geometry->Z0=" << geometry->z0 << std::endl;

  if(getVeryVerbose())
  {
    std::cout << ss.str() << std::endl;
  }
  // Now we actually initialize the data to something. If a subclass is involved
  // then the subclasses version of initializeData() will be used instead
  initializeData();

  //Doing a check sum to verify with matlab
  for (uint16_t y = 0; y < geometry->N_y; y++)
  {
    sum = 0;
    for (uint16_t x = 0; x < geometry->N_x; x++)
    {
      for (uint16_t z = 0; z < geometry->N_z; z++)
      {
//        sum += geometry->Object->d[k][j][i];
        sum += geometry->Object->getValue(z, x, y);
      }
    }
    ss << "Geometry check sum Y:" << y << " Value:" << sum << std::endl;
  }
  //End of check sum
  if (getVeryVerbose())
  {
    std::cout << ss.str() << std::endl;
  }
  setErrorCondition(0);
  setErrorMessage("");
  ss.str("");
  ss << getNameOfClass() << " Complete";
  notify(ss.str(), 0, UpdateProgressMessage);
}
