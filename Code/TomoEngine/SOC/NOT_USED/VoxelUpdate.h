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

#error

#ifndef VOXELUPDATE_H_
#define VOXELUPDATE_H_


#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/allocate.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/Common/randlib.h"
#include "TomoEngine/Filters/TomoFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"


typedef struct {
    int16_t Iter;
    uint16_t* cost_counter;
    Int32ArrayType::Pointer Counter;
    UInt8ImageType::Pointer VisitCount;
    AMatrixCol*** TempCol;
    AMatrixCol* VoxelLineResponse;
    ScaleOffsetParams* NuisanceParams;
    RealVolumeType::Pointer ErrorSino;
    RealVolumeType::Pointer Weight;
    UInt8ImageType::Pointer Mask;
    RealArrayType::Pointer cost;
    RNGVars* RandomNumber;
 //   uint8_t BOUNDARYFLAG[3][3][3];
    DATA_TYPE FILTER[3][3][3];
    DATA_TYPE THETA1;
    DATA_TYPE THETA2;
 //   DATA_TYPE NEIGHBORHOOD[3][3][3];
    DATA_TYPE V;
    DATA_TYPE MRF_P;
    DATA_TYPE SIGMA_X_P;
    FILE* Fp2;

} VoxelUpdateValues;
/*
 *
 */
class TomoEngine_EXPORT VoxelUpdate : public TomoFilter
{
  public:
    MXA_SHARED_POINTERS(VoxelUpdate)
    MXA_STATIC_NEW_MACRO(VoxelUpdate);
    MXA_STATIC_NEW_SUPERCLASS(TomoFilter, VoxelUpdate);
    MXA_TYPE_MACRO_SUPER(VoxelUpdate, TomoFilter)

    virtual ~VoxelUpdate();

    MXA_INSTANCE_PROPERTY(VoxelUpdateValues*, VoxelUpdateValues);

    virtual void execute();


  protected:
    VoxelUpdate();

    /**
     * @brief Finds the min and max of the neighborhood. This is required prior to calling solve()
     */
    void minMax(DATA_TYPE *low,DATA_TYPE *high);

    /**
     * @brief
     * @param ErrorSino
     * @param Weight
     * @return
     */
    DATA_TYPE computeCost(RealVolumeType::Pointer ErrorSino, RealVolumeType::Pointer Weight);

  private:

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


    VoxelUpdate(const VoxelUpdate&); // Copy Constructor Not Implemented
    void operator=(const VoxelUpdate&); // Operator '=' Not Implemented
};


#endif /* VOXELUPDATE_H_ */
