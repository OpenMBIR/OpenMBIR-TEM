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



#ifndef COMPUTATIONENGINE_H_
#define COMPUTATIONENGINE_H_

//-- C Headers
#include <stdio.h>

//-- MXA Headers
#include "MXA/Common/MXASetGetMacros.h"

#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/AbstractFilter.h"
#include "MBIRLib/Common/Observer.h"
#include "MBIRLib/Common/VoxelUpdateList.h"
#include "MBIRLib/GenericFilters/CostData.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/BrightField/BFQGGMRFPriorModel.h"

#include "MBIRLib/BrightField/BFForwardModel.h"
#include "MBIRLib/BrightField/BFAMatrixCol.h"

#include "MBIRLib/Common/EIMTime.h"
#define START_TIMER uint64_t startm = EIMTOMO_getMilliSeconds();
#define STOP_TIMER uint64_t stopm = EIMTOMO_getMilliSeconds();
#define PRINT_TIME(msg)\
  std::cout << indent << msg << ": " << ((double)stopm-startm)/1000.0 << " seconds" << std::endl;

#if 0
struct ImgIdx {
    int32_t xidx;
    int32_t zidx;
};

struct List
{
    int32_t NumElts;
    struct ImgIdx* Array;
};

void printList(struct List InList);
void maxList(struct List InList);
void minList(struct List InList);
#endif

/**
 * @class ReconstructionEngine ReconstructionEngine.h TomoEngine/SOC/ReconstructionEngine.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @author Singanallur Venkatakrishnan (Purdue University)
 * @version 1.0
 */
class MBIRLib_EXPORT ReconstructionEngine : public AbstractFilter
{

  public:
    MXA_SHARED_POINTERS(ReconstructionEngine);
    MXA_TYPE_MACRO(ReconstructionEngine);
    MXA_STATIC_NEW_MACRO(ReconstructionEngine);

    MXA_INSTANCE_PROPERTY(TomoInputsPtr, TomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, Sinogram)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)
    MXA_INSTANCE_PROPERTY(AdvancedParametersPtr, AdvParams)

    MXA_INSTANCE_PROPERTY(BFForwardModel::Pointer, ForwardModel);

    static void InitializeTomoInputs(TomoInputsPtr);
    static void InitializeSinogram(SinogramPtr);
    static void InitializeGeometry(GeometryPtr);
    static void InitializeAdvancedParams(AdvancedParametersPtr v);

    virtual ~ReconstructionEngine();


    /**
     * @brief overload from super class
     * @return
     */
    void execute();

    Real_t absMaxArray(std::vector<Real_t> &Array);


  protected:
    // Protect this constructor because we want to force the use of the other
    ReconstructionEngine();

    /**
     * @brief
     * @param errorSinogram
     * @param Weight
     * @return
     */
    int calculateCost(CostData::Pointer cost,
                      SinogramPtr sinogram,
                      GeometryPtr geometry,
                      RealVolumeType::Pointer ErrorSino,
                      BFQGGMRF::BFQGGMRF_Values* BFQGGMRF_Values);


    Real_t computeCost(SinogramPtr sinogram,
                       GeometryPtr geometry,
                       RealVolumeType::Pointer errorSinogram,
                       BFQGGMRF::BFQGGMRF_Values* BFQGGMRF_Values);

    //Updating voxels
    uint8_t updateVoxels(
        int16_t OuterIter,
        int16_t Iter,
        std::vector<BFAMatrixCol::Pointer> &TempCol,
        RealVolumeType::Pointer ErrorSino,
        std::vector<BFAMatrixCol::Pointer> &VoxelLineResponse,
        CostData::Pointer cost,
        BFQGGMRF::BFQGGMRF_Values* BFQGGMRF_Values,
        RealImageType::Pointer magUpdateMap,
        RealImageType::Pointer filtMagUpdateMap,
        UInt8Image_t::Pointer magUpdateMask,
        UInt8Image_t::Pointer m_VisitCount,
        Real_t PrevMagSum,
        uint32_t EffIterCount);

    void initializeROIMask(UInt8Image_t::Pointer Mask);


    Real_t roiVolumeSum(UInt8Image_t::Pointer Mask);

    uint8_t stopCriteria(RealImageType::Pointer magUpdateMap,UInt8Image_t::Pointer magUpdateMask, Real_t PrevMagSum, uint32_t EffIterCount);

    /**
     * Code to take the magnitude map and filter it with a hamming window
     * Returns the filtered magnitude map
     */
    void ComputeVSC(RealImageType::Pointer magUpdateMap,
                    RealImageType::Pointer filtMagUpdateMap);

    //Sort the entries of filtMagUpdateMap and set the threshold to be ? percentile
    Real_t SetNonHomThreshold(RealImageType::Pointer magUpdateMap);

    //Generate a regular sequential order list
    void GenRegularList(VoxelUpdateList::Pointer InpList);

    //Generating a list of indices for ICD updates. Output stored in the
    //private member m_VoxelIdxList
//    VoxelUpdateList::Pointer GenRandList(struct List InpList);

    //Generating a list for NHICD based on the magnitude update map
    VoxelUpdateList::Pointer GenNonHomList(Real_t NHThresh, RealImageType::Pointer magUpdateMap);

    /**
     * @brief
     */
    void calculateArithmeticMean();

    /**
     * @brief This is to be implemented at some point
     */
    void updateVoxelValues_NHICD();

    /**
     * @brief
     */
    int readInputData();

    /**
     * @brief
     */
    int initializeRoughReconstructionData();

    /**
     * @brief Calculates the x boundaries so the writers do not write extra data
     * @param cropStart
     * @param cropEnd
     */
    void computeOriginalXDims(uint16_t &cropStart, uint16_t &cropEnd);

    /**
     * @brief
     * @param H_t
     * @param OffsetT
     */
    void initializeHt(RealVolumeType::Pointer H_t, Real_t OffsetT);

    /**
     * @brief
     * @param H_t
     * @param VoxelLineResponse
     * @param haadfParameters
     */
    void storeVoxelResponse(RealVolumeType::Pointer H_t,
                            std::vector<BFAMatrixCol::Pointer> &VoxelLineResponse,
                            BFDetectorParameters::Pointer haadfParameters);

    /**
     * @brief
     * @param Y_Est
     * @param value
     */
    void initializeVolume(RealVolumeType::Pointer Y_Est, double value);

    /**
     * @brief
     * @param filepath
     */
    void writeReconstructionFile(const std::string &filepath);

    /**
     * @brief
     * @param vtkFile
     * @param cropStart
     * @param cropEnd
     */
    void writeVtkFile(const std::string &vtkFile, uint16_t cropStart, uint16_t cropEnd);

    /**
     * @brief
     * @param mrcFile
     * @param cropStart
     * @param cropEnd
     */
    void writeMRCFile(const std::string &mrcFile, uint16_t cropStart, uint16_t cropEnd);

    /**
     * @brief
     * @param file
     * @param cropStart
     * @param cropEnd
     */
    void writeAvizoFile(const std::string &file, uint16_t cropStart, uint16_t cropEnd);


    uint32_t Partition(RealArrayType::Pointer A,uint32_t p,uint32_t r);
    Real_t RandomizedSelect(RealArrayType::Pointer A,uint32_t p, uint32_t r,uint32_t i);
    uint32_t RandomizedPartition(RealArrayType::Pointer A,uint32_t p,uint32_t r);

  private:

    int m_NumThreads;

    VoxelUpdateList::Pointer m_VoxelIdxList;

    Real_t m_HammingWindow[5][5];

    /**
     * @brief
     */
    RealImageType::Pointer calculateVoxelProfile();

    /**
     * @brief
     * @return
     */
#ifndef EIMTOMO_USE_BFQGGMRF
    double surrogateFunctionBasedMin(Real_t currentVoxelValue);
#endif

    /**
    * @brief Updates a single line of voxels along y-axis
    * @param j_new
    * @param k_new
    */
    void updateVoxelLine(uint16_t j_new,uint16_t k_new);


    template<typename T>
    double solve(T* f, double a, double b, double err, int32_t *code,uint32_t iteration_count)

    {

      int signa, signb, signc;
      double fa, fb, fc, c; //, signaling_nan();
      double dist;
      uint32_t num_iter=0;


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

      while (num_iter < iteration_count)//(dist > err)
      {
        num_iter++;
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

    ReconstructionEngine(const ReconstructionEngine&); // Copy Constructor Not Implemented
    void operator=(const ReconstructionEngine&); // Operator '=' Not Implemented
};

#endif //CompEngine
