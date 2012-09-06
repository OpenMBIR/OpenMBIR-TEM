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

  // BrightField Data
    MXA_INSTANCE_PROPERTY(bool, UseBrightFieldData)
    MXA_INSTANCE_PROPERTY(TomoInputsPtr, BFTomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, BFSinogram)

    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialGain)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialOffset)
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, InitialVariance)

    MXA_INSTANCE_PROPERTY(Real_t, TargetGain)
    MXA_INSTANCE_PROPERTY(Real_t, DefaultOffset)
    MXA_INSTANCE_PROPERTY(Real_t, DefaultVariance)
    MXA_INSTANCE_PROPERTY(bool, UseDefaultOffset)

      // Bright Field Data Flag
    MXA_INSTANCE_PROPERTY(bool, BF_Flag)

      // These are the Nuisance Parameters that we need to solve for
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, I_0) //Gains
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Mu) //Offset
    MXA_INSTANCE_PROPERTY(RealArrayType::Pointer, Alpha) //Noise variance refinement factor

    MXA_INSTANCE_PROPERTY(RealVolumeType::Pointer, Weight) //This contains weights for each measurement = The diagonal covariance matrix in the Cost Func formulation

    void setQGGMRFValues(QGGMRF::QGGMRF_Values* qggmrf_values);

    void printNuisanceParameters(SinogramPtr m_Sinogram);

    /**
     * @brief Sets all the Nuisance Parameters Pointers to NULL.
     */
    void resetNuisanceParameters();

    void allocateNuisanceParameters(SinogramPtr m_Sinogram);

    void gainAndOffsetInitialization(uint16_t N_theta);

    void weightInitialization(size_t dims[3]);

    void costInitialization(SinogramPtr m_Sinogram);

  //  void initializePriorModel(TomoInputsPtr m_TomoInputs);

    int initializeBrightFieldData(SinogramPtr m_Sinogram);

    void initializeROIMask(SinogramPtr m_Sinogram, GeometryPtr m_Geometry, UInt8Image_t::Pointer Mask);

    int createNuisanceParameters(SinogramPtr m_Sinogram);

    int createInitialGainsData(SinogramPtr m_Sinogram);
    int createInitialOffsetsData(SinogramPtr m_Sinogram);
    int createInitialVariancesData(SinogramPtr m_Sinogram);

    /**
     *
     * @return
     */
    virtual int forwardProject(SinogramPtr m_Sinogram, GeometryPtr m_Geometry,
                               std::vector<HAADFAMatrixCol::Pointer> &TempCol,
                               std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
                               RealVolumeType::Pointer Y_Est,
                               RealVolumeType::Pointer ErrorSino);


    /**
     *
     */
    uint8_t updateVoxels(SinogramPtr m_Sinogram, GeometryPtr m_Geometry,
                         int16_t OuterIter, int16_t Iter,
                            UInt8Image_t::Pointer VisitCount,
                            std::vector<HAADFAMatrixCol::Pointer> &TempCol,
                            RealVolumeType::Pointer ErrorSino,
                            std::vector<HAADFAMatrixCol::Pointer> &VoxelLineResponse,
                            CostData::Pointer cost );

    /**
     *
     * @param m_Sinogram
     * @param ErrorSino
     * @param Y_Est
     */
    void calculateMeasurementWeight(SinogramPtr m_Sinogram, RealVolumeType::Pointer ErrorSino, RealVolumeType::Pointer Y_Est);
    /**
     *
     * @param m_Sinogram
     * @param ErrorSino
     * @param Y_Est
     * @param cost
     * @return
     */
    int jointEstimation(SinogramPtr m_Sinogram,
                        RealVolumeType::Pointer ErrorSino,
                        RealVolumeType::Pointer Y_Est, CostData::Pointer cost);
    /**
     *
     * @param cost
     * @param m_Sinogram
     * @param m_Geometry
     * @param ErrorSino
     * @return
     */
    int calculateCost(CostData::Pointer cost,
                      SinogramPtr m_Sinogram,
                      GeometryPtr m_Geometry,
                      RealVolumeType::Pointer ErrorSino,
                      QGGMRF::QGGMRF_Values *qggmrf_Values);
    /**
     * @brief
     * @param ErrorSino
     * @param Weight
     * @return
     */
    Real_t computeCost(SinogramPtr m_Sinogram,
                       GeometryPtr m_Geometry,
                       RealVolumeType::Pointer ErrorSino,
                       QGGMRF::QGGMRF_Values* qggmrf_Values);

    void updateWeights(SinogramPtr m_Sinogram,
                       RealVolumeType::Pointer ErrorSino);

    /**
     * Code to take the magnitude map and filter it with a hamming window
     * Returns the filtered magnitude map
     */
    void ComputeVSC(RealImageType::Pointer MagUpdateMap,
                    RealImageType::Pointer FiltMagUpdateMap,
                    GeometryPtr m_Geometry);

    //Sort the entries of FiltMagUpdateMap and set the threshold to be ? percentile
    Real_t SetNonHomThreshold(GeometryPtr m_Geometry, RealImageType::Pointer MagUpdateMap);


    void writeNuisanceParameters(SinogramPtr m_Sinogram);

    void writeSinogramFile(SinogramPtr m_Sinogram,
                           RealVolumeType::Pointer Final_Sinogram);

  protected:
    HAADFForwardModel();


  private:
    RealImageType::Pointer QuadraticParameters; //holds the coefficients of N_theta quadratic equations. This will be initialized inside the MAPICDREconstruct function
    RealImageType::Pointer Qk_cost;
    RealImageType::Pointer bk_cost;
    RealArrayType::Pointer ck_cost; //these are the terms of the quadratic cost function
    RealArrayType::Pointer d1;
    RealArrayType::Pointer d2; //hold the intermediate values needed to compute optimal mu_k


    Real_t HAMMING_WINDOW[5][5];
    //Markov Random Field Prior parameters - Globals DATA_TYPE
    Real_t FILTER[27];


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
