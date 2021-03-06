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



#ifndef _HAADF_ReconstructionEngine_H_
#define _HAADF_ReconstructionEngine_H_

//-- C Headers
#include <stdio.h>

//-- MXA Headers
#include "MXA/Common/MXASetGetMacros.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/AbstractFilter.h"
#include "MBIRLib/Common/Observer.h"
#include "MBIRLib/Common/AMatrixCol.h"
#include "MBIRLib/GenericFilters/CostData.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/HAADF/HAADFConstants.h"
#include "MBIRLib/HAADF/HAADF_ForwardModel.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/Reconstruction/QGGMRF_Functions.h"


/**
 * @class HAADF_ReconstructionEngine HAADF_ReconstructionEngine.h TomoEngine/SOC/HAADF_ReconstructionEngine.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @author Singanallur Venkatakrishnan (Purdue University)
 * @version 1.0
 */
class MBIRLib_EXPORT HAADF_ReconstructionEngine : public AbstractFilter
{

  public:
    MXA_SHARED_POINTERS(HAADF_ReconstructionEngine)
    MXA_TYPE_MACRO(HAADF_ReconstructionEngine)
    MXA_STATIC_NEW_MACRO(HAADF_ReconstructionEngine)

    MXA_INSTANCE_PROPERTY(TomoInputsPtr, TomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, Sinogram)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)
    MXA_INSTANCE_PROPERTY(AdvancedParametersPtr, AdvParams)

    MXA_INSTANCE_PROPERTY(HAADF_ForwardModel::Pointer, ForwardModel)

    MXA_INSTANCE_PROPERTY(bool, UseBrightFieldData)
    MXA_INSTANCE_PROPERTY(TomoInputsPtr, BFTomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, BFSinogram)

    static void InitializeTomoInputs(TomoInputsPtr);
    static void InitializeSinogram(SinogramPtr);
    static void InitializeGeometry(GeometryPtr);
    static void InitializeAdvancedParams(AdvancedParametersPtr v);

    static void InitializeScaleOffsetParams(HAADF_ForwardModel* forwardModel);

    virtual ~HAADF_ReconstructionEngine();

    /**
     * @brief overload from super class
     * @return
     */
    void execute();

    Real_t absMaxArray(std::vector<Real_t>& Array);

  protected:
    // Protect this constructor because we want to force the use of the other
    HAADF_ReconstructionEngine();

    void initVariables();

    int calculateCost(CostData::Pointer cost,
                      RealVolumeType::Pointer Weight,
                      RealVolumeType::Pointer ErrorSino);

    Real_t computeCost(RealVolumeType::Pointer ErrorSino,
                       RealVolumeType::Pointer Weight);


    /**
     *
     */
    uint8_t updateVoxels(int16_t OuterIter, int16_t Iter,
                         unsigned int updateType,
                         UInt8Image_t::Pointer VisitCount,
                         std::vector<AMatrixCol::Pointer>& TempCol,
                         RealVolumeType::Pointer ErrorSino,
                         RealVolumeType::Pointer Weight,
                         std::vector<AMatrixCol::Pointer>& VoxelLineResponse,
                         HAADF_ForwardModel* forwardModel,
                         UInt8Image_t::Pointer Mask,
                         CostData::Pointer cost);

    void initializeROIMask(UInt8Image_t::Pointer Mask);


    /**
     * Code to take the magnitude map and filter it with a hamming window
     * Returns the filtered magnitude map
     */
    void ComputeVSC();

    //Sort the entries of FiltMagUpdateMap and set the threshold to be ? percentile
    Real_t SetNonHomThreshold();


    void calculateGeometricMeanConstraint();

    /**
     * @brief
     */
    void calculateArithmeticMean();

    /**
     * @brief This is to be implemented at some point
     */
    void updateVoxelValues_NHICD();

    /**
     * @brief
     */
    int readInputData();

    /**
     * @brief
     */
    int initializeRoughReconstructionData();

    /**
     * @brief Calculates the x boundaries so the writers do not write extra data
     */
    void computeOriginalXDims(uint16_t& cropStart, uint16_t& cropEnd);

    /**
     * @brief
     * @param H_t
     * @param OffsetT
     */
//    void initializeHt(RealVolumeType::Pointer H_t, Real_t OffsetT);
//    void initializeHt(RealVolumeType::Pointer H_t);

    /**
     * @brief
     * @param H_t
     * @param VoxelLineResponse
     */
    void storeVoxelResponse(RealVolumeType::Pointer H_t,
                            std::vector<AMatrixCol::Pointer>& VoxelLineResponse);


    /**
     * @brief
     * @param Y_Est
     * @param value
     */
    void initializeVolume(RealVolumeType::Pointer Y_Est, double value);

    int initializeBrightFieldData();


    void gainAndOffsetInitialization();




    void calculateMeasurementWeight(RealVolumeType::Pointer Weight,
                                    RealVolumeType::Pointer ErrorSino,
                                    RealVolumeType::Pointer Y_Est);

    int jointEstimation(RealVolumeType::Pointer Weight,
                        RealVolumeType::Pointer ErrorSino,
                        RealVolumeType::Pointer Y_Est,
                        CostData::Pointer cost);

    void costInitialization(SinogramPtr sinogram);

    void updateWeights(RealVolumeType::Pointer Weight,
                       RealVolumeType::Pointer ErrorSino);

    int createNuisanceParameters(SinogramPtr sinogram);
    void printNuisanceParameters(SinogramPtr sinogram);
    void allocateNuisanceParameters();

#ifdef BF_RECON
    void processRawCounts();
#endif

  private:
    int m_NumThreads;

    //if 1 then this is NOT outside the support region; If 0 then that pixel should not be considered
    uint8_t BOUNDARYFLAG[27];
    //Markov Random Field Prior parameters - Globals DATA_TYPE
    Real_t FILTER[27];
    Real_t HAMMING_WINDOW[5][5];
    Real_t THETA1;
    Real_t THETA2;
    Real_t NEIGHBORHOOD[27];


#ifdef EIMTOMO_USE_QGGMRF
    QGGMRF::QGGMRF_Values m_QGGMRF_Values;
#else
    Real_t MRF_P;
    Real_t SIGMA_X_P;
#endif

    DetectorParameters::Pointer m_DetectorParameters;

    RealImageType::Pointer QuadraticParameters; //holds the coefficients of N_theta quadratic equations. This will be initialized inside the MAPICDREconstruct function

    RealImageType::Pointer MagUpdateMap;//Hold the magnitude of the reconstuction along each voxel line
    RealImageType::Pointer FiltMagUpdateMap;//Filters the above to compute threshold
    UInt8Image_t::Pointer MagUpdateMask;//Masks only the voxels of interest

    RealImageType::Pointer Qk_cost;
    RealImageType::Pointer bk_cost;
    RealArrayType::Pointer ck_cost; //these are the terms of the quadratic cost function
    RealArrayType::Pointer d1;
    RealArrayType::Pointer d2; //hold the intermediate values needed to compute optimal mu_k
    uint16_t NumOfViews; //this is kind of redundant but in order to avoid repeatedly send this info to the rooting function we save number of views
    Real_t LogGain; //again these information  are available but to prevent repeatedly sending it to the rooting functions we store it in a variable

    uint64_t startm;
    uint64_t stopm;




    /**
     * @brief
     * @param DetectorResponse
     * @param H_t
     * @return
     */
    RealVolumeType::Pointer forwardProject(RealVolumeType::Pointer DetectorResponse,
                                           RealVolumeType::Pointer H_t);

    /**
     * @brief
     * @param low
     * @param high
     */
    void minMax(Real_t* low, Real_t* high, Real_t currentVoxelValue);

    /**
     * @brief
     */
    RealImageType::Pointer calculateVoxelProfile();

    /**
     * @brief
     * @param row
     * @param col
     * @param slice
     * @param VoxelProfile
     */
    // void* calculateAMatrixColumn(uint16_t row, uint16_t col, uint16_t slice, DATA_TYPE** VoxelProfile);

    /**
     * @brief
     * @param ErrorSino
     * @param Weight
     * @return
     */
    Real_t estimateSigmaX(RealVolumeType::Pointer ErrorSino,
                          RealVolumeType::Pointer Weight);
    /**
     * @brief
     * @param row
     * @param col
     * @param slice
     * @param DetectorResponse
     */
    AMatrixCol::Pointer calculateAMatrixColumnPartial(uint16_t row, uint16_t col, uint16_t slice,
                                                      RealVolumeType::Pointer DetectorResponse);

    /**
     * @brief
     * @return
     */
#ifndef EIMTOMO_USE_QGGMRF
    double surrogateFunctionBasedMin(Real_t currentVoxelValue);
#endif


    /**
    * @brief Updates a single line of voxels along y-axis
    * @param j_new
    * @param k_new
    */
    void UpdateVoxelLine(uint16_t j_new, uint16_t k_new);


    template<typename T>
    double solve(T* f, double a, double b, double err, int32_t* code, uint32_t iteration_count)

    {

      int signa, signb, signc;
      double fa, fb, fc, c; //, signaling_nan();
      double dist;
      uint32_t num_iter = 0;


      fa = f->execute(a);
      signa = fa > 0;
      fb = f->execute(b);
      signb = fb > 0;

      /* check starting conditions */
      if(signa == signb)
      {
        if(signa == 1) { *code = 1; }
        else { *code = -1; }
        return (0.0);
      }
      else { *code = 0; }

      /* half interval search */
      if((dist = b - a) < 0) { dist = -dist; }

      while (num_iter < iteration_count)//(dist > err)
      {
        num_iter++;
        c = (b + a) / 2;
        fc = f->execute(c);
        signc = fc > 0;
        if(signa == signc)
        {
          a = c;
          fa = fc;
        }
        else
        {
          b = c;
          fb = fc;
        }
        if((dist = b - a) < 0) { dist = -dist; }
      }

      /* linear interpolation */
      if((fb - fa) == 0) { return (a); }
      else
      {
        c = (a * fb - b * fa) / (fb - fa);
        return (c);
      }
    }

    HAADF_ReconstructionEngine(const HAADF_ReconstructionEngine&); // Copy Constructor Not Implemented
    void operator=(const HAADF_ReconstructionEngine&); // Operator '=' Not Implemented
};

#endif //CompEngine
