/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Singanallur Venkatakrishnan (Purdue University)
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

// q-GGMRF is of the form
// \rho(\delta)= (|\deta/sigma|^p)/(c + |\deta/sigma|^(p-q))
// p=2; 1<=q<=2

#ifndef QGGMRF_FUNCTIONS_H_
#define QGGMRF_FUNCTIONS_H_

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"

namespace QGGMRF {


  static const unsigned int QGGMRF_ITER = 1;

  typedef struct {
      Real_t MRF_P;
      Real_t SIGMA_X_P;
      Real_t MRF_Q;
      Real_t MRF_C;
      Real_t MRF_ALPHA;
      Real_t SIGMA_X_P_Q;
      Real_t SIGMA_X_Q;
    Real_t gamma;
  } QGGMRF_Values;

  /**
   * @brief initializePriorModel
   * @param tomoInputs
   * @param qggmrf_values
   */
  void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF_Values* qggmrf_values);

/**
 *
 * @param delta
 * @param qggmrf_values
 * @return
 */
Real_t Value(Real_t delta, QGGMRF_Values* qggmrf_values);

/**
 *
 * @param delta
 * @param qggmrf_values
 * @return
 */
Real_t Derivative(Real_t delta, QGGMRF_Values* qggmrf_values);

/**
 *
 * @param delta
 * @param qggmrf_values
 * @return
 */
Real_t SecondDerivative(Real_t delta, QGGMRF_Values* qggmrf_values);

/**
 *
 * @param umin
 * @param umax
 * @param RefValue
 * @param BOUNDARYFLAG
 * @param NEIGHBORHOOD
 * @param qggmrf_values
 * @param QGGMRF_Params
 */
void ComputeParameters(Real_t umin, Real_t umax, Real_t RefValue,
                       uint8_t* BOUNDARYFLAG, Real_t* NEIGHBORHOOD,
                       QGGMRF::QGGMRF_Values* qggmrf_values,
                       Real_t* QGGMRF_Params);
/**
 *
 * @param umin
 * @param umax
 * @param currentVoxelValue
 * @param BOUNDARYFLAG
 * @param FILTER
 * @param NEIGHBORHOOD
 * @param THETA1
 * @param THETA2
 * @param qggmrf_values
 * @return
 */
Real_t FunctionalSubstitution(Real_t umin, Real_t umax, Real_t currentVoxelValue,
                                 uint8_t* BOUNDARYFLAG, Real_t* FILTER, Real_t* NEIGHBORHOOD,
                                 Real_t THETA1, Real_t THETA2,
                                 QGGMRF_Values* qggmrf_values);



Real_t updatePriorModel(Real_t NewSigmaX,QGGMRF_Values* qggmrf_values);
} /* end namespace QGGMRF */

#endif /* QGGMRF_FUNCTIONS_H_ */
