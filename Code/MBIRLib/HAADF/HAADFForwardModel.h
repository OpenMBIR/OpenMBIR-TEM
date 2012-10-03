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

/*
 * The gains,offsets, variance (the very existence of those) + calculation of the
 * AMatrixColumn + Voxel Line Response ~= Forward Model specification
 *
 *
 * So are the Weights. The way we set the weight is a function of the forward model
 * too. And the nuisance parameters are also part of the forward model.
 *
 * The forward model is used in a bunch of ways. Initially it forward projects the initial
 * volume and computes a error between the data and the "forward projection".

 * Also in the iterations it is used to compute the various parameters (THETA1 and THETA2)
 * both of which require scanning through a portion of the the AMatrix, Voxel Line response,
 * Nuisance parameters , weights, bright field data etc. Which reminds me the BrightField
 * (which we have disable in the GUI is also a part of the forward model if its available)
 *
 * Not to the AMatrix and Voxel Line but to the nuisance parameters. So in some way
 * our forward model has some "unknown" things hence forcing us to update some values
 * of it. It like having an equation some coefficients I know , some I am going to figure out.
 *
 * So "Forward Project" is the algorithm and the "Forward Model" is the input into
 * the algorithm and the output from the algorithm is a reconstructed volume that
 * we compare to the original MRC data which gives us the "Error Sinogram"?
 *
 * Right. The output is not a volume but a tilt series (so we can compare it to
 * the acquired data).
 */

#ifndef _HAADFFORWARDMODEL_H_
#define _HAADFFORWARDMODEL_H_

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/AbstractFilter.h"
#include "MBIRLib/Common/Observer.h"
#include "MBIRLib/IOFilters/MRCWriter.h"

#include "MXA/Common/MXASetGetMacros.h"
#include "MBIRLib/Common/Observable.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"

#include "MBIRLib/GenericFilters/CostData.h"
#include "MBIRLib/HAADF/QGGMRFPriorModel.h"
#include "MBIRLib/HAADF/HAADFAMatrixCol.h"


/**
 * @class HAADFForwardModel HAADFForwardModel.h MBIRLib/Reconstruction/HAADFForwardModel.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Sep 5, 2012
 * @version 1.0
 */
class HAADFForwardModel : public Observable
{
  public:

    MXA_SHARED_POINTERS(HAADFForwardModel)
    MXA_TYPE_MACRO(HAADFForwardModel)
    MXA_STATIC_NEW_MACRO(HAADFForwardModel)

    virtual ~HAADFForwardModel();

    MXA_INSTANCE_PROPERTY(bool, Verbose)
    MXA_INSTANCE_PROPERTY(bool, VeryVerbose)
    MXA_INSTANCE_PROPERTY(int, ErrorCondition)
    MXA_INSTANCE_PROPERTY(bool, Cancel)

    MXA_INSTANCE_PROPERTY(AdvancedParametersPtr, AdvParams)
    MXA_INSTANCE_PROPERTY(TomoInputsPtr, TomoInputs)

    // Bright Field Data Flag
    MXA_INSTANCE_PROPERTY(bool, BF_Flag)
    MXA_INSTANCE_PROPERTY(TomoInputsPtr, BFTomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, BFSinogram)

    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialGain)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialOffset)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialVariance)

    MXA_INSTANCE_PROPERTY(Real_t, TargetGain)
    MXA_INSTANCE_PROPERTY(Real_t, DefaultOffset)
    MXA_INSTANCE_PROPERTY(Real_t, DefaultVariance)
    MXA_INSTANCE_PROPERTY(bool, UseDefaultOffset)
	MXA_INSTANCE_PROPERTY(Real_t, BraggThreshold)
    


      // These are the Nuisance Parameters that we need to solve for
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, I_0) //Gains
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Mu) //Offset
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Alpha) //Noise variance refinement factor

    MXA_INSTANCE_PROPERTY(RealVolumeType::Pointer, Weight) //This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation
    MXA_INSTANCE_PROPERTY(UInt8VolumeType::Pointer, Selector); //This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation
	
    void setQGGMRFValues(QGGMRF::QGGMRF_Values* qggmrf_values);

    void printNuisanceParameters(SinogramPtr sinogram);

    /**
     * @brief Sets all the Nuisance Parameters Pointers to NULL.
     */
    void resetNuisanceParameters();

    void allocateNuisanceParameters(SinogramPtr sinogram);

    void gainAndOffsetInitialization(uint16_t N_theta);

    void weightInitialization(size_t dims[3]);

    void costInitialization(SinogramPtr sinogram);

    int initializeBrightFieldData(SinogramPtr sinogram);

    int createNuisanceParameters(SinogramPtr sinogram);

    int createInitialGainsData(SinogramPtr sinogram);
    int createInitialOffsetsData(SinogramPtr sinogram);
    int createInitialVariancesData(SinogramPtr sinogram);

    /**
     *
     * @return
     */
    virtual int forwardProject(SinogramPtr sinogram, GeometryPtr geometry,
                               std::vector<HAADFAMatrixCol::Pointer> &TempCol,
                               std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
                               RealVolumeType::Pointer yEstimate,
                               RealVolumeType::Pointer errorSinogram);


    /**
     *
     */
  

    /**
     *
     * @param sinogram
     * @param errorSinogram
     * @param yEstimate
     */
    void calculateMeasurementWeight(SinogramPtr sinogram, RealVolumeType::Pointer errorSinogram, RealVolumeType::Pointer yEstimate);
    /**
     *
     * @param sinogram
     * @param errorSinogram
     * @param yEstimate
     * @param cost
     * @return
     */
    void jointEstimation(SinogramPtr sinogram,
                        RealVolumeType::Pointer errorSinogram,
                        RealVolumeType::Pointer yEstimate, CostData::Pointer cost);
	void updateWeights(SinogramPtr sinogram,
                       RealVolumeType::Pointer errorSinogram);
#ifdef BF_RECON
	void updateSelector(SinogramPtr sinogram,
                       RealVolumeType::Pointer errorSinogram);
	
	void printRatioSelected(SinogramPtr sinogram);
	
	void writeSelectorMrc(SinogramPtr sinogram,GeometryPtr geometry,RealVolumeType::Pointer ErrorSino);
	
	Real_t estimateBraggThresold(SinogramPtr sinogram,RealVolumeType::Pointer ErrorSino,Real_t percentage);
	//void setUpBraggThreshold(Real_t Threshold);
	
#endif

  
    void writeNuisanceParameters(SinogramPtr sinogram);

    void writeSinogramFile(SinogramPtr sinogram,
                           RealVolumeType::Pointer finalSinogram);
	
	Real_t forwardCost(SinogramPtr sinogram,
					   RealVolumeType::Pointer errorSinogram);
	
	//Computing the theta parameters of the cost function
	void computeTheta(size_t Index,
						 std::vector<HAADFAMatrixCol::Pointer> &TempCol,
						 int32_t xzSliceIdx,
						 std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
						 RealVolumeType::Pointer ErrorSino,
						 SinogramPtr sinogram,
					     Int32ArrayType::Pointer Thetas);
	
	//After updating a voxel "Index",update sinogram
	void updateErrorSinogram(Real_t ChangeInVoxelValue,
							 size_t Index,
							 std::vector<HAADFAMatrixCol::Pointer> &TempCol,
							 int32_t xzSliceIdx,
							 std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
							 RealVolumeType::Pointer errorSinogram,
							 SinogramPtr sinogram);
	
#ifdef BF_RECON
	void processRawCounts(SinogramPtr sinogram);
#endif

  protected:
    HAADFForwardModel();

  private:
    RealImageType::Pointer m_QuadraticParameters; //holds the coefficients of N_theta quadratic equations. This will be initialized inside the MAPICDREconstruct function
    RealImageType::Pointer m_QkCost;
    RealImageType::Pointer m_BkCost;
    RealArrayType::Pointer m_CkCost; //these are the terms of the quadratic cost function
    RealArrayType::Pointer m_D1;
    RealArrayType::Pointer m_D2; //hold the intermediate values needed to compute optimal mu_k

    Real_t k_HammingWindow[5][5];
	Real_t Theta[2]; //Theta1 and Theta2 in the optimization
 
#ifdef EIMTOMO_USE_QGGMRF
    QGGMRF::QGGMRF_Values* m_QGGMRF_Values;
#else
    Real_t MRF_P;
    Real_t SIGMA_X_P;
#endif


    HAADFForwardModel(const HAADFForwardModel&); // Copy Constructor Not Implemented
    void operator=(const HAADFForwardModel&); // Operator '=' Not Implemented
};

#endif /* _HAADFFORWARDMODEL_H_ */
