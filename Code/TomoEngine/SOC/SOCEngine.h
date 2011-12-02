/*
 *  ComputationEngineOptimized.h
 *  HAADFSTEM_Xcode
 *
 *  Created by Singanallur Venkatakrishnan on 6/29/11.
 *  Copyright 2011 Purdue University. All rights reserved.
 *
 */

#ifndef COMPUTATIONENGINE_H_
#define COMPUTATIONENGINE_H_

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractPipeline.h"
#include "TomoEngine/Common/Observer.h"

#include "TomoEngine/SOC/SOCStructures.h"

/**
 *
 */
class SOCEngine : public AbstractPipeline, public Observer
{

  public:
    MXA_SHARED_POINTERS(SOCEngine);
    MXA_TYPE_MACRO(SOCEngine);
    MXA_STATIC_NEW_MACRO(SOCEngine);

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);


    virtual ~SOCEngine();

    /**
     * @brief overload from super class
     * @return
     */
    void execute();



    DATA_TYPE absMaxArray(std::vector<DATA_TYPE> &Array, uint16_t NumElts);



  protected:
    // Protect this constructor because we want to force the use of the other
    SOCEngine();

    void initVariables();

    void initializeSinoParameters();

    void initializeGeomParameters();


  private:

    uint8_t BOUNDARYFLAG[3][3][3]; //if 1 then this is NOT outside the support region; If 0 then that pixel should not be considered
    //Markov Random Field Prior parameters - Globals DATA_TYPE
    DATA_TYPE FILTER[3][3][3];
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

    DATA_TYPE* cosine;
    DATA_TYPE* sine; //used to store cosine and sine of all angles through which sample is tilted
    DATA_TYPE* BeamProfile; //used to store the shape of the e-beam
    DATA_TYPE BEAM_WIDTH;
    DATA_TYPE OffsetR;
    DATA_TYPE OffsetT;

    DATA_TYPE** QuadraticParameters; //holds the coefficients of N_theta quadratic equations. This will be initialized inside the MAPICDREconstruct function
    DATA_TYPE** Qk_cost;
    DATA_TYPE** bk_cost;
    DATA_TYPE* ck_cost; //these are the terms of the quadratic cost function
    DATA_TYPE* d1;
    DATA_TYPE* d2; //hold the intermediate values needed to compute optimal mu_k
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
    DATA_TYPE*** forwardProject(DATA_TYPE*** DetectorResponse, DATA_TYPE*** H_t);

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
    void* calculateVoxelProfile();

    /**
     * @brief
     * @param row
     * @param col
     * @param slice
     * @param VoxelProfile
     */
    void* calculateAMatrixColumn(uint16_t row, uint16_t col, uint16_t slice, DATA_TYPE** VoxelProfile);

    /**
     * @brief
     * @param ErrorSino
     * @param Weight
     * @return
     */
    DATA_TYPE computeCost(DATA_TYPE*** ErrorSino, DATA_TYPE*** Weight);

    /**
     * @brief
     * @param row
     * @param col
     * @param slice
     * @param DetectorResponse
     */
    void* calculateAMatrixColumnPartial(uint16_t row, uint16_t col, uint16_t slice, DATA_TYPE*** DetectorResponse);

    /**
     * @brief
     * @param row
     * @param col
     * @param VoxelProfile
     */
    void* detectorResponse(uint16_t row, uint16_t col, DATA_TYPE** VoxelProfile);

    /**
     * @brief
     * @return
     */
    double surrogateFunctionBasedMin();

    /**
     * @brief  Solves equation (*f)(x) = 0 on x in [a,b]. Uses half interval method.
     * Requires that (*f)(a) and (*f)(b) have opposite signs.
     * Returns code=0 if signs are opposite.
     * Returns code=1 if signs are both positive.
     * Returns code=-1 if signs are both negative.
     * @param f  pointer to object that implements the function to be solved
     * @param a minimum value of solution
     * @param b maximum value of solution
     * @param err accuarcy of solution
     * @param code error code
     * @return
     */
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
