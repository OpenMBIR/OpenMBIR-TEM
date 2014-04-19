/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
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
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef TOMOARRAY_HPP_
#define TOMOARRAY_HPP_

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "MBIRLib/Reconstruction/ReconstructionConstants.h"

/**
 * @brief Creates a new Array by allocating memory in a contiguous space
 * @param dims The dimensions of your data with the SLOWEST moving dimension
 * in the first (zero) slot and continuing to the FASTEST moving dimension.
 * @param name The human name for the array. This can be helpful for debugging
 * @return Shared Pointer to the data
 *
 *
 * Example: If you have a 2D array that you would like to allocate and the
 * dimensions are a width (X) of 512 pixels and a height (Y) of 256 pixels
 * and you want the data laid out so that when traversing the array one would
 * walk the width of the image FIRST then the heigth then you would create
 * a TomoArray2D like the following:
 * @code
 *  size_t dims[2] = { 256, 512 };
 *  TomoArray2D<int, int*, 2>::Pointer image = TomoArray2D<int, int*, 2>::New(dims, "Image);
 *
 *    Width (X) --->
 *    - - - - - - - - -
 * H | | | | | | | | | |
 * e  - - - - - - - - -
 * i | | | | | | | | | |
 * g  - - - - - - - - -
 * h | | | | | | | | | |
 * t  - - - - - - - - -
 *   | | | | | | | | | |
 * Y  - - - - - - - - -

 * @endcode
 * @n
 * Example: If you have a 3D array that you would like to allocate and the
 * dimensions are a width (X) of 512 pixels, height (Y) of 256 pixels and depth (Z)
 * of 128 pixels and you want the data laid out so that when traversing the
 * array one would walk the width of the image FIRST, then the heigth, then
 * the Depth then you would create a TomoArray3D like the following:
 * @code
 *  size_t dims[2] = { 128, 256, 512 };
 *  TomoArray3D<int, int*, 3>::Pointer image = TomoArray3D<int, int*, 3>::New(dims, "Image);
 *
 *    Width (X) --->
 *    - - - - - - - - -
 * H | | | | | | | | | |
 * e  - - - - - - - - -
 * i | | | | | | | | | |
 * g  - - - - - - - - -
 * h | | | | | | | | | |
 * t  - - - - - - - - -
 *   | | | | | | | | | |
 * Y  - - - - - - - - -

 * @endcode
 */


template<typename T, typename Ptr, int SIZE>
class TomoArray
{
  public:
    typedef TomoArray<T, Ptr, SIZE> Self;
    typedef boost::shared_ptr<Self> Pointer;
    typedef boost::shared_ptr<const Self> ConstPointer;
    typedef boost::weak_ptr<TomoArray<T, Ptr, SIZE> > WeakPointer;
    typedef boost::weak_ptr<TomoArray<T, Ptr, SIZE> > ConstWeakPointer;

    /** This is the raw pointer to the data.
     */
    Ptr d;

    /**
     * @brief Creates a NULL wrapped pointer
     * @return
     */
    static Pointer NullPointer(void)
    {
      return Pointer(static_cast<TomoArray<T, Ptr, SIZE>*>(0));
    }

    /**
     * @brief Creates a new Array by allocating memory in a contiguous space
     * @param dims The dimensions of your data with the SLOWEST moving dimension
     * in the first (zero) slot and continuing to the FASTEST moving dimension.
     * @param name The human name for the array. This can be helpful for debugging
     * @return Shared Pointer to the data
     */
    static Pointer New(size_t* dims, const std::string& name)
    {
      //assert(SIZE < 4);

      Pointer sharedPtr(new TomoArray<T, Ptr, SIZE>(dims));
      sharedPtr->setName(name);
      return sharedPtr;
    }

    virtual ~TomoArray()
    {
#if 0
      if (m_Name.empty() == false)
      {
        std::cout << "Deallocating TomoArray " << m_Name << std::endl;
      }
#endif
      if (SIZE == 1 || SIZE == 2 || SIZE == 3)
      {
        free(d);
      }
      else
      {
        deallocate((T*)d, SIZE);
      }
    }

    inline size_t numElements()
    {
      size_t count = 1;
      for(int i = 0; i < SIZE; ++i)
      {
        count += m_Dims[i];
      }
      return count;
    }

    inline void initializeWithZeros()
    {
      size_t count = 1;
      for(int i = 0; i < SIZE; ++i)
      {
        count += m_Dims[i];
      }
      ::memset(d, 0, count * sizeof(T));
    }

    /* ******************* These are 3D array methods ********************* */
    inline size_t calcIndex(size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      return (m_Dims[1] * m_Dims[2] * z) + (m_Dims[2] * y) + (x);
    }

    inline T getValue(size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      return d[(m_Dims[1] * m_Dims[2] * z) + (m_Dims[2] * y) + (x)];
    }

    inline void setValue(T v, size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      d[(m_Dims[1]*m_Dims[2]*z) + (m_Dims[2]*y) + (x)] = v;
    }

    inline void addToValue(T v, size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      d[(m_Dims[1]*m_Dims[2]*z) + (m_Dims[2]*y) + (x)] += v;
    }

    inline void deleteFromValue(T v, size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      d[(m_Dims[1]*m_Dims[2]*z) + (m_Dims[2]*y) + (x)] -= v;
    }

    inline void divideByValue(T v, size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      d[(m_Dims[1]*m_Dims[2]*z) + (m_Dims[2]*y) + (x)] /= v;
    }

    inline void multiplyByValue(T v, size_t z, size_t y, size_t x)
    {
      assert(SIZE == 3);
      d[(m_Dims[1]*m_Dims[2]*z) + (m_Dims[2]*y) + (x)] *= v;
    }

    inline size_t calcIndex(size_t slow, size_t fast)
    {
      assert(SIZE == 2);
      return (m_Dims[1] * slow) + (fast);
    }

    /* ******************* These are 2D array methods ********************* */
    inline T getValue(size_t slow, size_t fast)
    {
      assert(SIZE == 2);
      return d[(m_Dims[1] * slow) + (fast)];
    }

    inline void setValue(T v, size_t slow, size_t fast)
    {
      assert(SIZE == 2);
      d[(m_Dims[1]*slow) + (fast)] = v;
    }

    inline T* getPointer(size_t slow, size_t fast)
    {
      assert(SIZE == 2);
      return d + (m_Dims[1] * slow) + (fast);
    }


    /* ******************* These are 1D array methods ********************* */
    inline T getValue(size_t x)
    {
      return d[x];
    }

    inline void setValue(T v, size_t x)
    {
      d[x] = v;
    }

    inline T* getPointer(size_t offset)
    {
      return d + offset;
    }


    void setName(const std::string& name) { m_Name = name;}
    Ptr getPointer() { return d; }
    size_t* getDims() {return m_Dims; }
    int getNDims() { return m_NDims; }
    int getTypeSize() { return sizeof(T); }

  protected:
    TomoArray(size_t* dims)
    {
      size_t total = 1;
      for(size_t i = 0; i < SIZE; ++i)
      {
        m_Dims[i] = dims[i];
        total *= m_Dims[i];
      }
      if (SIZE == 1 || SIZE == 2 || SIZE == 3)
      {
        d = reinterpret_cast<Ptr>(malloc(sizeof(T) * total));
      }
      else
      {
        d = reinterpret_cast<Ptr>(allocate(sizeof(T), SIZE, m_Dims));
      }
      m_NDims = SIZE;
    }

  private:
    size_t m_Dims[SIZE];
    int  m_NDims;
    std::string m_Name;


    T* allocate(size_t s, unsigned int d, size_t* d1)
    {
      // va_list ap;             /* varargs list traverser */
      size_t max,                /* size of array to be declared */
             *q;                     /* pointer to dimension list */
      char** r,               /* pointer to beginning of the array of the
                                   * pointers for a dimension */
           **s1, *t, *tree;        /* base pointer to beginning of first array */
      size_t i, j;               /* loop counters */
      assert(false);
      // int *d1;                /* dimension list */

      // va_start(ap,d);
      // d1 = (int *) mget_spc(d,sizeof(int));

      //      for(i=0;i<d;i++)
      //        d1[i] = va_arg(ap,int);

      r = &tree;
      q = d1;                /* first dimension */
      max = 1;
      for (i = 0; i < d - 1; i++, q++)
      {
        /* for each of the dimensions
                                                 * but the last */
        max *= (*q);
        //        r[0]=(char *)mget_spc(max,sizeof(char**));
        r[0] = (char*)malloc((size_t)(max * sizeof(char**)));
        //   printf("max: %d   sizeof(char **): %d  r: %p r[0]: 0x%08x \n", max , sizeof(char **), r, r[0]);
        r = (char**) r[0];     /* step through to beginning of next
                                 * dimension array */
      }
      max *= s * (*q);        /* grab actual array memory */
      //r[0] = (char*)mget_spc(max,sizeof(char));
      r[0] = (char*)malloc((size_t)(max * sizeof(char)));
      // printf("max: %d   sizeof(char): %d  r: %p r[0]: 0x%08x \n", max , sizeof(char), r, r[0]);
      /*
       * r is now set to point to the beginning of each array so that we can
       * use it to scan down each array rather than having to go across and
       * then down
       */
      r = (char**) tree;      /* back to the beginning of list of arrays */
      q = d1;                 /* back to the first dimension */
      max = 1;
      for (i = 0; i < d - 2; i++, q++)
      {
        /* we deal with the last
                                                 * array of pointers later on */
        max *= (*q);    /* number of elements in this dimension */
        for (j = 1, s1 = r + 1, t = r[0]; j < max; j++)
        {
          /* scans down array for
                                                   * first and subsequent
                                                   * elements */

          /*  modify each of the pointers so that it points to
          * the correct position (sub-array) of the next
          * dimension array. s1 is the current position in the
          * current array. t is the current position in the
          * next array. t is incremented before s1 is, but it
          * starts off one behind. *(q+1) is the dimension of
          * the next array. */

          *s1 = (t += sizeof (char**)** (q + 1));
          s1++;
        }
        r = (char**) r[0];     /* step through to begining of next
                                 * dimension array */
      }
      max *= (*q);              /* max is total number of elements in the
                                 * last pointer array */

      /* same as previous loop, but different size factor */
      for (j = 1, s1 = r + 1, t = r[0]; j < max; j++)
      { *s1++ = (t += s** (q + 1)); }

      //      va_end(ap);
      //      free((void *)d1);
      return (T*)tree;              /* return base pointer */
    }
    /*
     * multifree releases all memory that we have already declared analogous to
     * free() when using malloc()
     */

    void deallocate(T* r, int d)
    {
      void** p;
      void* next = NULL;
      int i;

      for (p = (void**)r, i = 0; i < d; p = (void**)next, i++)
      {
        if(p != NULL)
        {
          next = *p;
          // printf("p: %p  Value of P: 0x%08x  i:  %d   Value of Next: 0x%08x \n", p, *p, i, next);
          free((void*)p);
        }
      }
    }

    TomoArray(const TomoArray&); // Copy Constructor Not Implemented
    void operator=(const TomoArray&); // Operator '=' Not Implemented

};

typedef TomoArray<uint8_t, uint8_t*, 3> UInt8VolumeType;
typedef TomoArray<uint8_t, uint8_t*, 2> UInt8Image_t;
typedef TomoArray<uint8_t, uint8_t*, 1> UInt8ArrayType;

typedef TomoArray<int32_t, int32_t*, 1> Int32ArrayType;
typedef TomoArray<uint32_t, uint32_t*, 1> UInt32ArrayType;

typedef TomoArray<Real_t, Real_t*, 3> RealVolumeType;
typedef TomoArray<Real_t, Real_t*, 2> RealImageType;
typedef TomoArray<Real_t, Real_t*, 1> RealArrayType;

typedef TomoArray<float, float*, 3> FloatVolumeType;
typedef TomoArray<float, float*, 2> FloatImageType;
typedef TomoArray<float, float*, 1> FloatArrayType;

#endif /* TOMOARRAY_HPP_ */
