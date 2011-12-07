/*
 * WrappedArrayTest.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include <stdio.h>
#include <iostream>


#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/TomoArray.hpp"



int main(int argc, char **argv)
{

  typedef TomoArray<int, 3>   Int3DArrayType;



  int x = 3;
  int y = 6;
  int z = 2;
  int dims[3] = { x, y, z};

  {
    Int3DArrayType::Pointer object = Int3DArrayType::New(dims);

    int*** data = reinterpret_cast<int***>(object->getPointer());
    printf("Array to Hold Pointers to start of Slice: %p \n", data);
    for(int i = 0; i < z; ++i)
    {
      printf("| 0x%08x ", data[i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    printf("Array to Hold Pointers to start of each Row: %p \n", &(data[0][0]));
    for(int i =0; i < y*z; ++i)
    {
      printf("| 0x%08x ", data[0][i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    for(int k = 0; k < z; ++k)
      {
        printf("Slice %d  Ptr: %p   Value: 0x%08x \n", k, data[k], *data[k]);
        for(int j = 0; j < y; ++j)
        {
          printf("%d %p | ", j, data[k][j]);
          for(int i=0; i < x; ++i)
          {
            data[k][j][i] = i*10 + j*100 + k*1000;
            printf("{[%d][%d][%d]  %04d %p}  ", k, j, i, data[k][j][i], &(data[k][j][i]));
          }
          printf("\n");
        }
      }
    printf("\n");
    std::cout << "Cleaning up Memory" << std::endl;
    object = Int3DArrayType::NullPointer();
  }

  typedef TomoArray<double, 2> Double2DArrayType;
  {
    {
      Double2DArrayType::Pointer object = Double2DArrayType::New(dims);

      double*** data = reinterpret_cast<double***>(object->getPointer());
      printf("Array to Hold Pointers to start of each Row: %p \n", data);
      for(int i = 0; i < y; ++i)
      {
        printf("| 0x%08x ", data[i]);
      }
      printf("|\n");
      printf("-------------------------- \n");

      printf("Data Array Address: %p Num. Elements: %d\n", &(data[0][0]), x * y);
      printf("-------------------------- \n");

      std::cout << "Cleaning up Memory" << std::endl;
      object = Double2DArrayType::NullPointer();
    }
  }



  return 0;
}


