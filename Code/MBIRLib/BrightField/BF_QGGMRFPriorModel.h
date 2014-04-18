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

//q-GGMRF is of the form
// \rho(\delta)= (|\deta/sigma|^p)/(c + |\deta/sigma|^(p-q))
// p=2; 1<=q<=2

#ifndef _BF_QGGMRFPriorModel_H_
#define _BF_QGGMRFPriorModel_H_

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Reconstruction/QGGMRF_Functions.h"

namespace QGGMRF {


  /**
   * @brief initializePriorModel
   * @param tomoInputs
   * @param qggmrf_values
   */
  void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF::QGGMRF_Values* qggmrf_values, Real_t* filter);

  /**
   * @brief initializeFilter
   * @param filter
   */
  void initializeFilter(Real_t* filter);

  /**
   * @brief PriorModelCost
   * @param geometry
   * @param qggmrf_values
   * @param filter
   * @return
   */
  Real_t PriorModelCost(GeometryPtr geometry, QGGMRF::QGGMRF_Values* qggmrf_values, Real_t* filter);


} /* end namespace BFQGGMRF */

#endif /* BFQGGMRF_FUNCTIONS_H_ */
