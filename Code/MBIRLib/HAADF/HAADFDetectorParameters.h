/*
 * HAADFParameters.h
 *
 *  Created on: Sep 5, 2012
 *      Author: mjackson
 */

#ifndef _HAADFPARAMETERS_H_
#define _HAADFPARAMETERS_H_

#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"

/*
 *
 */
class HAADFDetectorParameters
{
  public:
    MXA_SHARED_POINTERS(HAADFDetectorParameters);
    MXA_STATIC_NEW_MACRO(HAADFDetectorParameters)
    virtual ~HAADFDetectorParameters();
    //used to store cosine and sine of all angles through which sample is tilted
    MXA_INSTANCE_PROPERTY( RealArrayType::Pointer, cosine)
    MXA_INSTANCE_PROPERTY( RealArrayType::Pointer, sine)
    MXA_INSTANCE_PROPERTY( RealArrayType::Pointer, BeamProfile) //used to store the shape of the e-beam
    MXA_INSTANCE_PROPERTY( Real_t, BEAM_WIDTH)
    MXA_INSTANCE_PROPERTY( Real_t, OffsetR)
    MXA_INSTANCE_PROPERTY( Real_t, OffsetT)

    void calculateSinCos(SinogramPtr m_Sinogram);
    void initializeBeamProfile(SinogramPtr m_Sinogram, AdvancedParametersPtr m_AdvParams);

  protected:
    HAADFDetectorParameters();


  private:




    HAADFDetectorParameters(const HAADFDetectorParameters&); // Copy Constructor Not Implemented
    void operator=(const HAADFDetectorParameters&); // Operator '=' Not Implemented
};

#endif /* _HAADFPARAMETERS_H_ */
