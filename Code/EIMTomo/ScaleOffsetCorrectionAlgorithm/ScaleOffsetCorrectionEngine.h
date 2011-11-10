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

#include "EIMTomo/EIMTomo.h"


#include "ScaleOffsetStructures.h"
#include "ScaleOffsetCorrectionInputs.h"


/**
 *
 */
class SOCEngine
{

  public:
    SOCEngine();
    virtual ~SOCEngine();

    DATA_TYPE*** ForwardProject(Sino*, Geom*, DATA_TYPE***, DATA_TYPE***);
    void CE_CalculateSinCos(Sino*);
    void CE_InitializeBeamProfile(Sino*);
    void CE_MinMax(DATA_TYPE*, DATA_TYPE*);
    void* CE_CalculateVoxelProfile(Sino*, Geom*);
    int CE_MAPICDReconstruct(Sino*, Geom*, CommandLineInputs*);
    void* CE_CalculateAMatrixColumn(uint16_t, uint16_t, uint16_t, Sino*, Geom*, DATA_TYPE**);

    double CE_DerivOfCostFunc(double);

    DATA_TYPE CE_ComputeCost(DATA_TYPE***, DATA_TYPE***, Sino*, Geom*);
    void* CE_CalculateAMatrixColumnPartial(uint16_t, uint16_t, uint16_t, Sino*, Geom*, DATA_TYPE***);
    void* CE_DetectorResponse(uint16_t, uint16_t, Sino*, Geom*, DATA_TYPE**);

    template<class functor>
    double solve(functor* f, /* pointer to function to be solved */
                 double a, /* minimum value of solution */
                 double b, /* maximum value of solution */
                 double err, /* accuarcy of solution */
                 int32_t *code /* error code */
                 )
    /* Solves equation (*f)(x) = 0 on x in [a,b]. Uses half interval method.*/
    /* Requires that (*f)(a) and (*f)(b) have opposite signs.   */
    /* Returns code=0 if signs are opposite.        */
    /* Returns code=1 if signs are both positive.         */
    /* Returns code=-1 if signs are both negative.        */
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

    double CE_SurrogateFunctionBasedMin();
    //  double CE_ConstraintEquation(DATA_TYPE);
    // DATA_TYPE* CE_RootsOfQuadraticFunction(DATA_TYPE, DATA_TYPE, DATA_TYPE);
#ifdef QGGMRF
    DATA_TYPE CE_QGGMRF_Value(DATA_TYPE );
    DATA_TYPE CE_QGGMRF_Derivative(DATA_TYPE);
    DATA_TYPE CE_QGGMRF_SecondDerivative(DATA_TYPE);
    void CE_ComputeQGGMRFParameters(DATA_TYPE ,DATA_TYPE);
    DATA_TYPE CE_FunctionalSubstitution(DATA_TYPE ,DATA_TYPE );
#endif//qggmrf
  private:

    bool CE_Cancel;
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

    SOCEngine(const SOCEngine&); // Copy Constructor Not Implemented
    void operator=(const SOCEngine&); // Operator '=' Not Implemented
};

#endif //CompEngine
