#ifndef _HAADF_QGGMRF_PRIORMODEL_H_
#define _HAADF_QGGMRF_PRIORMODEL_H_

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/HAADF/HAADFConstants.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"

class QGGMRF_Values;

namespace QGGMRF {


  /**
   * @brief initializePriorModel
   * @param tomoInputs
   * @param qggmrf_values
   */
  void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF_Values* qggmrf_values);

  /**
   * @brief updatePriorModel
   * @param NewSigmaX
   * @param qggmrf_values
   * @return
   */
  Real_t updatePriorModel(Real_t NewSigmaX, QGGMRF_Values* qggmrf_values, Real_t gamma);
}

#endif /* _HAADF_QGGMRF_PRIORMODEL_H_ */
