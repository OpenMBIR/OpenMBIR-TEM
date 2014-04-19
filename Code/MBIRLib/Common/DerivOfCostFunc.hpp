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
#ifndef DERIVOFCOSTFUNC_HPP_
#define DERIVOFCOSTFUNC_HPP_



/**
 * @brief
 * @author
 * @version
 */
class DerivOfCostFunc
{

  public:
    DerivOfCostFunc(uint8_t _BOUNDARYFLAG[3][3][3],
                    Real_t _NEIGHBORHOOD[3][3][3],
                    Real_t _FILTER[3][3][3],
                    Real_t _V,
                    Real_t _THETA1,
                    Real_t _THETA2,
                    Real_t _SIGMA_X_P,
                    Real_t _MRF_P) :

      V(_V),
      THETA1(_THETA1),
      THETA2(_THETA2),
      SIGMA_X_P(_SIGMA_X_P),
      MRF_P(_MRF_P)
    {

      for(size_t i = 0; i < 3; ++i)
      {
        for(size_t j = 0; j < 3; ++j)
        {
          for(size_t k = 0; k < 3; ++k)
          {
            BOUNDARYFLAG[i][j][k] = _BOUNDARYFLAG[i][j][k];
            NEIGHBORHOOD[i][j][k] = _NEIGHBORHOOD[i][j][k];
            FILTER[i][j][k] = _FILTER[i][j][k];
          }
        }
      }
    }

    virtual ~DerivOfCostFunc()
    {
    }

    double execute(Real_t u) const
    {
      double temp = 0;
      double value = 0;
      uint8_t i, j, k;
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          for (k = 0; k < 3; k++)
          {
            if(BOUNDARYFLAG[i][j][k] == 1)
            {
              if(u - NEIGHBORHOOD[i][j][k] >= 0.0)
              { temp += ((double)FILTER[i][j][k] * (1.0) * pow(fabs(u - (double)NEIGHBORHOOD[i][j][k]), (double)(MRF_P - 1))); }
              else { temp += ((double)FILTER[i][j][k] * (-1.0) * pow(fabs(u - (double)NEIGHBORHOOD[i][j][k]), (double)(MRF_P - 1))); }
            }
          }

      //printf("V While updating %lf\n",V);
      //scanf("Enter value %d\n",&k);
      value = THETA1 + THETA2 * (u - V) + (temp / SIGMA_X_P);

      return value;
    }

  protected:
    DerivOfCostFunc() {}

  private:
    uint8_t BOUNDARYFLAG[3][3][3];
    Real_t NEIGHBORHOOD[3][3][3];
    Real_t FILTER[3][3][3];
    Real_t V;
    Real_t THETA1;
    Real_t THETA2;
    Real_t SIGMA_X_P;
    Real_t MRF_P;

    DerivOfCostFunc(const DerivOfCostFunc&); // Copy Constructor Not Implemented
    void operator=(const DerivOfCostFunc&); // Operator '=' Not Implemented
};


#endif /* DERIVOFCOSTFUNC_HPP_ */
