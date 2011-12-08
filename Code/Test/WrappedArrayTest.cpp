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
#include "TomoEngine/SOC/SOCStructures.h"


int main(int argc, char **argv)
{





  size_t x = 3;
  size_t y = 6;
  size_t z = 2;
  size_t dims[3] = { x, y, z};
  std::cout << "X, Y, Z: " << x << " " << y << " " << z << std::endl;
  std::cout << "-----------------Testing 3D Array------------------" << std::endl;
  {
    Int32VolumeType::Pointer object = Int32VolumeType::New(dims);

    int32_t*** data = object->d;

    printf("Array to Hold Pointers to start of Slice: Address=%p  Num Elements:%d \n Values:", data, z);
    for(int i = 0; i < z; ++i)
    {
      printf("| 0x%08x ", object->d[i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    printf("Array to Hold Pointers to start of each Row: Address=%p  Num Elements:%d \n Values:", &(data[0][0]), z*y);
    for(int i =0; i < y*z; ++i)
    {
      printf("| 0x%08x ", object->d[0][i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    for(int k = 0; k < z; ++k)
      {
        printf("Slice %d  Ptr: %p   Value: 0x%08x \n", k, object->d[k], *object->d[k]);
        for(int j = 0; j < y; ++j)
        {
          printf("Row %d %p [%d][%d][0] | ", j, object->d[k][j], k, j);
          for(int i=0; i < x; ++i)
          {
            object->d[k][j][i] = i*10 + j*100 + k*1000;
            printf(" %04d |", data[k][j][i]);
          }
          printf("\n");
        }
      }
    printf("\n");
    std::cout << "Cleaning up Memory\n" << std::endl;
    object = Int32VolumeType::NullPointer();
  }


  std::cout << "-----------------Testing 2D Array------------------" << std::endl;

  {

    RealImageType::Pointer object = RealImageType::New(dims);

    DATA_TYPE** data = reinterpret_cast<double**>(object->getPointer());
    printf("Array to Hold Pointers to start of each Row: %p \n", data);
    for(int i = 0; i < y; ++i)
    {
      printf("| 0x%08x ", data[i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    printf("Data Array Address: %p Num. Elements: %d\n", &(data[0][0]), x * y);
    printf("-------------------------- \n");

    std::cout << "Cleaning up Memory\n" << std::endl;
    object = RealImageType::NullPointer();

  }

  std::cout << "-----------------Testing 1D Array------------------" << std::endl;
  {
    RealArrayType::Pointer object = RealArrayType::New(dims);
    DATA_TYPE* data = object->getPointer();
    printf("Array to Hold Pointers to start of Data: Address=%p Num.Elements: %d\n", data, x);

  }


  return 0;
}


