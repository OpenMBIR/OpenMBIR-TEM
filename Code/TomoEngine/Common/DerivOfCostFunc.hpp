/*
 * DerivOfCostFunc.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: mjackson
 */

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
                    DATA_TYPE _NEIGHBORHOOD[3][3][3],
                    DATA_TYPE _FILTER[3][3][3],
                    DATA_TYPE _V,
                    DATA_TYPE _THETA1,
                    DATA_TYPE _THETA2,
                    DATA_TYPE _SIGMA_X_P,
                    DATA_TYPE _MRF_P) :

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

    double execute(DATA_TYPE u) const
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
          temp +=((double)FILTER[i][j][k] * (1.0) * pow(fabs(u - (double)NEIGHBORHOOD[i][j][k]), (double)(MRF_P - 1)));
              else temp += ((double)FILTER[i][j][k] * (-1.0) * pow(fabs(u - (double)NEIGHBORHOOD[i][j][k]), (double)(MRF_P - 1)));
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
    DATA_TYPE NEIGHBORHOOD[3][3][3];
    DATA_TYPE FILTER[3][3][3];
    DATA_TYPE V;
    DATA_TYPE THETA1;
    DATA_TYPE THETA2;
    DATA_TYPE SIGMA_X_P;
    DATA_TYPE MRF_P;

    DerivOfCostFunc(const DerivOfCostFunc&); // Copy Constructor Not Implemented
    void operator=(const DerivOfCostFunc&); // Operator '=' Not Implemented
};


#endif /* DERIVOFCOSTFUNC_HPP_ */
