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

#ifndef CE_CONSTRAINTEQUATION_HPP_
#define CE_CONSTRAINTEQUATION_HPP_

#include <limits>

#include "TomoEngine/TomoEngine.h"

/**
 * @class CE_ConstraintEquation CE_ConstraintEquation.h EIMTomo/ScaleOffsetCorrectionEngine.h
 * @brief
 * @author
 * @date Nov 1, 2011
 * @version 1.0
 */
class CE_ConstraintEquation
{
  public:
    CE_ConstraintEquation(uint16_t NumOfViews,
                          RealImageType::Pointer QuadraticParameters,
                          RealArrayType::Pointer d1,
                          RealArrayType::Pointer d2,
                          RealImageType::Pointer Qk_cost,
                          RealImageType::Pointer bk_cost,
                          RealArrayType::Pointer ck_cost,
                          DATA_TYPE LogGain) :
        NumOfViews(NumOfViews),
        QuadraticParameters(QuadraticParameters),
        d1(d1),
        d2(d2),
        Qk_cost(Qk_cost),
        bk_cost(bk_cost),
        ck_cost(ck_cost),
        LogGain(LogGain)
    {
    }//constructor which assigns values to the private members of the class ; The members of the class and the parameters have the same name here

    virtual ~CE_ConstraintEquation()
    {
    }

    /**
     *
     * @param a
     * @param b
     * @param c
     * @return
     */
    DATA_TYPE* CE_RootsOfQuadraticFunction(DATA_TYPE a, DATA_TYPE b, DATA_TYPE c)
    {
      DATA_TYPE* temp = (DATA_TYPE*)get_spc(2, sizeof(double));
      DATA_TYPE value = 0, discriminant;
      temp[0] = -1;
      temp[1] = -1;
      discriminant = b * b - 4 * a * c;
      if(discriminant < 0)
      {
        printf("Quadratic has no real roots\n");
        return temp;
      }
      else
      {
        value = sqrt(discriminant);
        temp[0] = (-b + value) / (2 * a);
        temp[1] = (-b - value) / (2 * a);
      }
      return temp;

    }

    /**
     *
     * @param lamda
     * @return
     */
    double execute(DATA_TYPE lambda)
    {
      double sum = 0, temp_cost = 0, min = std::numeric_limits<double>::max();
      double value = 0;
      DATA_TYPE* root;
      double temp_mu;
      uint8_t i, min_index = 0;
      uint16_t i_theta;

      for (i_theta = 0; i_theta < NumOfViews; i_theta++)
      {
    //    root = CE_RootsOfQuadraticFunction(QuadraticParameters->getValue(i_theta,0), QuadraticParameters->getValue(i_theta,1), lambda);
        root = CE_RootsOfQuadraticFunction(QuadraticParameters->d[i_theta][0], QuadraticParameters->d[i_theta][1], lambda);
        //Evaluate which root results in a lower cost function
        for (i = 0; i < 2; i++)
        {
          if(root[i] > 0) // If the value of I0[k] is positive
          {
            temp_mu = d1->d[i_theta] - root[i] * d2->d[i_theta]; //for a given lambda we can calculate I0(\lambda) and hence mu(lambda)

//            temp_cost = (Qk_cost->getValue(i_theta, 0) * root[i] * root[i]
//                        +  2 * Qk_cost->getValue(i_theta, 1) * root[i] * temp_mu
//                        +  temp_mu * temp_mu * Qk_cost->getValue(i_theta, 2)
//                        - 2 * (bk_cost->getValue(i_theta, 0) * root[i] + temp_mu * bk_cost->getValue(i_theta, 1))
//                        + ck_cost->d[i_theta]); //evaluating the cost function


            temp_cost = (Qk_cost->d[i_theta][0] * root[i] * root[i]
                        +  2 * Qk_cost->d[i_theta][1] * root[i] * temp_mu
                        +  temp_mu * temp_mu * Qk_cost->d[i_theta][2]
                        - 2 * (bk_cost->d[i_theta][0] * root[i] + temp_mu * bk_cost->d[i_theta][1])
                        + ck_cost->d[i_theta]); //evaluating the cost function
            if(temp_cost < min)
            {
              min = temp_cost;
              min_index = i;

            }
          }
        }

        if(root[min_index] > 0) sum += log(root[min_index]); //max{(-b+sqrt(b^2 - 4*a*c))/2*a,(-b+sqrt(b^2 - 4*a*c))/2*a}
        else
        {
          printf("Log of a negative number\n");
          printf("View %d\n", i_theta);
          printf("Roots of the quadratic are %lf %lf \n", root[0], root[1]);
        }
        free(root);

      }
      value = sum - LogGain;
      return value;
    }

  protected:
    CE_ConstraintEquation() {}

  private:
    uint16_t NumOfViews;
    RealImageType::Pointer  QuadraticParameters;
    RealArrayType::Pointer d1;
    RealArrayType::Pointer d2; //hold the intermediate values needed to compute optimal mu_k
    RealImageType::Pointer Qk_cost;
    RealImageType::Pointer bk_cost;
    RealArrayType::Pointer ck_cost; //these are the terms of the quadratic cost function
    DATA_TYPE LogGain;
    CE_ConstraintEquation(const CE_ConstraintEquation&); // Copy Constructor Not Implemented
    void operator=(const CE_ConstraintEquation&); // Operator '=' Not Implemented
};



#endif /* CE_CONSTRAINTEQUATION_HPP_ */
