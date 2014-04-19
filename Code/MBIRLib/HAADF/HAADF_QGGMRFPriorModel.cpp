

#include "HAADF_QGGMRFPriorModel.h"
#include "MBIRLib/Reconstruction/QGGMRF_Functions.h"


namespace QGGMRF
{

  void initializePriorModel(TomoInputsPtr tomoInputs, QGGMRF_Values* qggmrf_values)
  {
    qggmrf_values->MRF_P = 2;
    qggmrf_values->MRF_Q = tomoInputs->p;
    qggmrf_values->MRF_C = 0.001;
    qggmrf_values->MRF_ALPHA = 1.5;
    qggmrf_values->SIGMA_X_P = pow(tomoInputs->SigmaX, qggmrf_values->MRF_P);
    qggmrf_values->SIGMA_X_P_Q = pow(tomoInputs->SigmaX, (qggmrf_values->MRF_P - qggmrf_values->MRF_Q));
    qggmrf_values->SIGMA_X_Q = pow(tomoInputs->SigmaX, qggmrf_values->MRF_Q);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  Real_t updatePriorModel(Real_t NewSigmaX, QGGMRF_Values* qggmrf_values, Real_t gamma)
  {
    Real_t SigmaX = gamma * pow(NewSigmaX * qggmrf_values->MRF_Q, 1.0 / qggmrf_values->MRF_Q);
    qggmrf_values->SIGMA_X_P = pow(SigmaX, qggmrf_values->MRF_P);
    qggmrf_values->SIGMA_X_P_Q = pow(SigmaX, (qggmrf_values->MRF_P - qggmrf_values->MRF_Q));
    qggmrf_values->SIGMA_X_Q = pow(SigmaX, qggmrf_values->MRF_Q);
    return SigmaX;
  }


} /* End Namespace */
