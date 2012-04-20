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
#include "TomoArray_Test.hpp"

void otherTest()
{
  size_t x = 3;
  size_t y = 5;
  size_t z = 2;
  size_t dims[3] = { x, y, z};
  std::cout << "X, Y, Z: " << x << " " << y << " " << z << std::endl;
  std::cout << "-----------------Testing 3D Array------------------" << std::endl;
  {
    Int32VolumeType::Pointer object = Int32VolumeType::New(dims, "Test3D");

    int32_t*** data = object->d;

    printf("Dim[0]=%d Dim[1]=%d Dim[2]=%d\n", object->dims[0], object->dims[1], object->dims[2]);

    printf("Array to Hold Pointers to start of Slice: Address=%p  Num Elements:%d \n Values:", data, dims[2]);
    for(size_t i = 0; i < dims[2]; ++i)
    {
      printf("| 0x%08x ", object->d[i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    printf("Array to Hold Pointers to start of each Row: Address=%p  Num Elements:%d \n Values:", &(data[0][0]), dims[2]*dims[1]);
    for(size_t i =0; i < dims[2]*dims[1]; ++i)
    {
      printf("| 0x%08x ", object->d[0][i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    ::memset( &(data[0][0][0]), 0xAB, sizeof(int32_t) * dims[0] * dims[1] * dims[2]);
    int32_t* ptr = &(data[0][0][0]);
    size_t idx = 0;
    for(size_t k = 0; k < z; ++k)
      {
        printf("Slice %lu  Ptr: %p   Value: 0x%08x \n", k, object->d[k], *object->d[k]);
        for(size_t j = 0; j < y; ++j)
        {
          printf("Row %lu %p [%lu][%lu][0] | ", j, object->d[k][j], k, j);
          for(size_t i=0; i < x; ++i)
          {
            object->d[k][j][i] = i*10 + j*100 + k*1000;
            printf(" %04d |", data[k][j][i]);
//            MAKE_3D_INDEX(idx, dims[1], dims[2], i,j,k);
//            assert(object->d[k][j][i] == ptr[idx]);
//            idx = object->calcIndex(k, j, i);
//            idx = object->calcIndex(i, j, k);
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
   // size_t idx = MAKE_2D_INDEX(QuadraticParameters->dims[1], i_theta, 0);

    RealImage_t::Pointer object = RealImage_t::New(dims, "Test2D");

    DATA_TYPE* data = reinterpret_cast<double*>(object->getPointer());
    printf("Array to Hold Pointers to start of each Row: %p \n", data);
    for(size_t i = 0; i < y; ++i)
    {
      printf("| 0x%08x ", data[i]);
    }
    printf("|\n");
    printf("-------------------------- \n");

    printf("Data Array Address: %p Num. Elements: %ld\n", data, x * y);
    printf("-------------------------- \n");

    std::cout << "Cleaning up Memory\n" << std::endl;
    object = RealImage_t::NullPointer();

  }

  std::cout << "-----------------Testing 1D Array------------------" << std::endl;
  {
    RealArrayType::Pointer object = RealArrayType::New(dims, "1D");
    DATA_TYPE* data = object->getPointer();
    printf("Array to Hold Pointers to start of Data: Address=%p Num.Elements: %lud\n", data, x);

  }

}


typedef TomoArray_Test<int32_t, int32_t**, 2> Int32ImageTest_t;

int main(int argc, char **argv)
{
  otherTest();



  return 0;
}


