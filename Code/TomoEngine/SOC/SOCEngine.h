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



#ifndef COMPUTATIONENGINE_H_
#define COMPUTATIONENGINE_H_

#include <stdio.h>

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractPipeline.h"
#include "TomoEngine/Common/Observer.h"
#include "TomoEngine/mt/mt19937ar.h"

#include "TomoEngine/SOC/SOCStructures.h"
#include "TomoEngine/Filters/CostData.h"


/**
 * @class SOCEngine SOCEngine.h TomoEngine/SOC/SOCEngine.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @author Singanallur Venkatakrishnan (Purdue University)
 * @version 1.0
 */
class TomoEngine_EXPORT SOCEngine : public AbstractPipeline, public Observer
{

  public:
    MXA_SHARED_POINTERS(SOCEngine);
    MXA_TYPE_MACRO(SOCEngine);
    MXA_STATIC_NEW_MACRO(SOCEngine);

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);
    MXA_INSTANCE_PROPERTY(ScaleOffsetParams*, NuisanceParams);


    virtual ~SOCEngine();

    enum VoxelUpdateType {
      RegularRandomOrderUpdate = 0,
      HomogeniousUpdate = 1,
      NonHomogeniousUpdate = 2
    };

    /**
     * @brief overload from super class
     * @return
     */
    void execute();

    DATA_TYPE absMaxArray(std::vector<DATA_TYPE> &Array);

  protected:
    // Protect this constructor because we want to force the use of the other
    SOCEngine();

    void initVariables();

    void calculateGeometricMeanConstraint(ScaleOffsetParams* NuisanceParams);

    /**
     * @brief This is to be implemented at some point
     */
    void updateVoxelValues_NHICD();

    void updateVoxels(int16_t Iter, VoxelUpdateType updateType, UInt8ImageType::Pointer VisitCount,
                      RNGVars* RandomNumber, AMatrixCol*** TempCol,
                      RealVolumeType::Pointer ErrorSino, RealVolumeType::Pointer Weight,
                      AMatrixCol* VoxelLineResponse, ScaleOffsetParams &NuisanceParams,
                      UInt8ImageType::Pointer Mask, CostData::Pointer cost);

  private:
    //if 1 then this is NOT outside the support region; If 0 then that pixel should not be considered
    uint8_t BOUNDARYFLAG[3][3][3];
    //Markov Random Field Prior parameters - Globals DATA_TYPE
    DATA_TYPE FILTER[3][3][3];
    DATA_TYPE HAMMING_WINDOW[5][5];
    DATA_TYPE THETA1;
    DATA_TYPE THETA2;
    DATA_TYPE NEIGHBORHOOD[3][3][3];
    DATA_TYPE V;
    DATA_TYPE MRF_P;
    DATA_TYPE SIGMA_X_P;
#ifdef QGGMRF
    //QGGMRF extras
    DATA_TYPE MRF_Q,MRF_C;
    DATA_TYPE QGGMRF_Params[26][3];
    DATA_TYPE MRF_ALPHA;
#endif
    //used to store cosine and sine of all angles through which sample is tilted
    RealArrayType::Pointer cosine;
    RealArrayType::Pointer sine;
    RealArrayType::Pointer BeamProfile; //used to store the shape of the e-beam
    DATA_TYPE BEAM_WIDTH;
    DATA_TYPE OffsetR;
    DATA_TYPE OffsetT;

    RealImageType::Pointer QuadraticParameters; //holds the coefficients of N_theta quadratic equations. This will be initialized inside the MAPICDREconstruct function

    RealImageType::Pointer MagUpdateMap;//Hold the magnitude of the reconstuction along each voxel line
    RealImageType::Pointer FiltMagUpdateMap;//Filters the above to compute threshold
    Uint8ImageType::Pointer MagUpdateMask;//Masks only the voxels of interest

    RealImageType::Pointer Qk_cost;
    RealImageType::Pointer bk_cost;
    RealArrayType::Pointer ck_cost; //these are the terms of the quadratic cost function
    RealArrayType::Pointer d1;
    RealArrayType::Pointer d2; //hold the intermediate values needed to compute optimal mu_k
    uint16_t NumOfViews; //this is kind of redundant but in order to avoid repeatedly send this info to the rooting function we save number of views
    DATA_TYPE LogGain; //again these information  are available but to prevent repeatedly sending it to the rooting functions we store it in a variable

    uint64_t startm;
    uint64_t stopm;

    /**
     * @brief
     * @param DetectorResponse
     * @param H_t
     * @return
     */
    RealVolumeType::Pointer forwardProject(RealVolumeType::Pointer DetectorResponse, RealVolumeType::Pointer H_t);

    /**
     * @brief
     */
    void calculateSinCos();

    /**
     * @brief
     */
    void initializeBeamProfile();

    /**
     * @brief
     * @param low
     * @param high
     */
    void minMax(DATA_TYPE *low, DATA_TYPE *high);

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
    DATA_TYPE computeCost(RealVolumeType::Pointer ErrorSino, RealVolumeType::Pointer Weight);

    /**
     * @brief
     * @param row
     * @param col
     * @param slice
     * @param DetectorResponse
     */
    void* calculateAMatrixColumnPartial(uint16_t row,uint16_t col, uint16_t slice, RealVolumeType::Pointer DetectorResponse);

    /**
     * @brief
     * @return
     */
    double surrogateFunctionBasedMin();

	//Updates a single line of voxels along y-axis
	void UpdateVoxelLine(uint16_t j_new,uint16_t k_new);


	/**
     * Code to take the magnitude map and filter it with a hamming window
     * Returns the filtered magnitude map
     */
    void ComputeVSC();

	//Sort the entries of FiltMagUpdateMap and set the threshold to be ? percentile
	DATA_TYPE SetNonHomThreshold();

    /**
     *Code to return the threshold corresponding to the top T percentage of magnitude
    */

    DATA_TYPE ComputeThreshold(RealImageType::Pointer FilteredMagnitudeMap);

    template<typename T>
    double solve(T* f, double a, double b, double err, int32_t *code)

    {
      int signa, signb, signc;
      double fa, fb, fc, c, signaling_nan();
      double dist;

      fa = f->execute(a);
      signa = fa > 0;
      fb = f->execute(b);
      signb = fb > 0;

      /* check starting conditions */
      if(signa == signb)
      {
        if(signa == 1) *code = 1;
        else *code = -1;
        return (0.0);
      }
      else *code = 0;

      /* half interval search */
      if((dist = b - a) < 0) dist = -dist;
      while (dist > err)
      {
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
        if((dist = b - a) < 0) dist = -dist;
      }

      /* linear interpolation */
      if((fb - fa) == 0) return (a);
      else
      {
        c = (a * fb - b * fa) / (fb - fa);
        return (c);
      }
    }

    SOCEngine(const SOCEngine&); // Copy Constructor Not Implemented
    void operator=(const SOCEngine&); // Operator '=' Not Implemented
};

#endif //CompEngine
