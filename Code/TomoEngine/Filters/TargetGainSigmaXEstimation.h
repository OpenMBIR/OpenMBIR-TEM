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



#ifndef TARGETGAINSIGMAXESTIMATION_H_
#define TARGETGAINSIGMAXESTIMATION_H_


#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


/**
 * @class TargetGainSigmaXEstimation TargetGainSigmaXEstimation.h TomoEngine/Filters/TargetGainSigmaXEstimation.h
 * @brief This filter will estimate the Target Gain and Sigma X parameters for a
 * given input MRC file.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Mar 13, 2012
 * @version 1.0
 */
class TomoEngine_EXPORT TargetGainSigmaXEstimation : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(TargetGainSigmaXEstimation)
    MXA_STATIC_NEW_MACRO(TargetGainSigmaXEstimation);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, TargetGainSigmaXEstimation);
    MXA_TYPE_MACRO_SUPER(TargetGainSigmaXEstimation, TomoFilter)

    virtual ~TargetGainSigmaXEstimation();

    MXA_INSTANCE_STRING_PROPERTY(InputFile);
    MXA_INSTANCE_PROPERTY(Real_t, SigmaXEstimate);
    MXA_INSTANCE_PROPERTY(Real_t, TargetGainEstimate);
    MXA_INSTANCE_PROPERTY(Real_t, SampleThickness);
    MXA_INSTANCE_PROPERTY(Real_t, DefaultOffset);
    MXA_INSTANCE_PROPERTY(unsigned int, TiltAngles);

    virtual void execute();

  protected:
    TargetGainSigmaXEstimation();

  private:


    template<typename T>
    void calcMinMax(T* data, int total, Real_t &min, Real_t &max, Real_t &sum2)
    {
      for (int i = 0; i < total; i++)
      {
          if(data[i] > max) max = data[i];
          if(data[i] < min) min = data[i];
          sum2 += (data[i] - m_DefaultOffset);
      }
    }


    TargetGainSigmaXEstimation(const TargetGainSigmaXEstimation&); // Copy Constructor Not Implemented
    void operator=(const TargetGainSigmaXEstimation&); // Operator '=' Not Implemented
};




#endif /* TARGETGAINSIGMAXESTIMATION_H_ */
