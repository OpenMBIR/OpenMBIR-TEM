/*
 * ScaleOffsetCorrectionConstants.h
 *
 *  Created on: Nov 11, 2011
 *      Author: mjackson
 */

#ifndef SCALEOFFSETCORRECTIONCONSTANTS_H_
#define SCALEOFFSETCORRECTIONCONSTANTS_H_

#include <string>



namespace ScaleOffsetCorrection
{

//  const std::string OutputDirectory("ScaleOffsetCorrection_Output");
  const std::string DetectorResponseFile("DetectorResponse.bin");
  const std::string CostFunctionFile("CostFunc.bin");
  const std::string FinalGainParametersFile("FinalGainParameters.bin");
  const std::string FinalOffsetParametersFile("FinalOffsetParameters.bin");
  const std::string FinalVariancesFile("FinalVariances.bin");
  const std::string ReconstructedSinogramFile("ReconstructedSino.bin");
  const std::string VoxelProfileFile("VoxelProfile.bin");
  const std::string ForwardProjectedObjectFile("ForwardProjectedObject.bin");
  const std::string CostFunctionCoefficientsFile("CostFunctionCoefficients.bin");

}



#endif /* SCALEOFFSETCORRECTIONCONSTANTS_H_ */
