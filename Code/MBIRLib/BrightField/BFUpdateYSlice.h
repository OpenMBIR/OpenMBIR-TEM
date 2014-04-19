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


#ifndef UPDATEYSLICE_H_
#define UPDATEYSLICE_H_


#include <vector>


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"
#include "MBIRLib/BrightField/BFConstants.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/BrightField/BFReconstructionEngine.h"
#include "MBIRLib/Common/AMatrixCol.h"
#include "MBIRLib/BrightField/BFForwardModel.h"



#define USE_TBB_TASK_GROUP 1
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#include <tbb/task.h>

#endif


//Updates a line of voxels
//Required Global vars / Private members
//ErrorSino
//Weights
//m_Geometry->Object - This is alrady available
//TempCol
//VoxelLine
//NuisanceParams


/* - These get initialized to Zero each time so they are essentially local variables */
//NEIGHBORHOOD
//BOUNDARYFLAG

/**
 * @brief Updates one or more Y Slices of voxels. This can be run serial or in a multi-threaded
 * fashion using the Threading Building Blocks library
 *
 */
class BFUpdateYSlice
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
  : public tbb::task
#endif
{
  public:
    /**
    * @brief
    * @param yStart
    * @param yEnd
    * @param geometry
    * @param outerIter
    * @param innerIter
    * @param sinogram
    * @param tempCol
    * @param errorSino
    * @param voxelLineResponse
    * @param forwardModel
    * @param mask
    * @param magUpdateMap
    * @param magUpdateMask
    * @param voxelUpdateType
    * @param averageUpdate
    * @param averageMagnitudeOfRecon
    * @param zeroSkipping
    * @param qggmrf_values
    */
    BFUpdateYSlice(uint16_t yStart, uint16_t yEnd,
                   GeometryPtr geometry, int16_t outerIter, int16_t innerIter,
                   SinogramPtr  sinogram,
                   std::vector<AMatrixCol::Pointer>& tempCol,
                   RealVolumeType::Pointer errorSino,
                   std::vector<AMatrixCol::Pointer>& voxelLineResponse,
                   BFForwardModel* forwardModel,
                   UInt8Image_t::Pointer mask,
                   RealImageType::Pointer magUpdateMap, //Hold the magnitude of the reconstuction along each voxel line
                   UInt8Image_t::Pointer magUpdateMask,
                   unsigned int voxelUpdateType,
                   Real_t* averageUpdate,
                   Real_t* averageMagnitudeOfRecon,
                   unsigned int zeroSkipping,
                   QGGMRF::QGGMRF_Values* qggmrf_values,
                   VoxelUpdateList::Pointer voxelUpdateList);

    ~BFUpdateYSlice();

    /**
     * @brief
     * @return
     */
    int getZeroCount();

    /**
    *
    */
#if defined (OpenMBIR_USE_PARALLEL_ALGORITHMS)
    tbb::task*
#else
    void
#endif
    execute();

  protected:

    /**
    * @brief
    */
    void initVariables();

  private:
    uint16_t m_YStart;
    uint16_t m_YEnd;
    GeometryPtr m_Geometry;
    int16_t m_OuterIter;
    int16_t m_InnerIter;
    SinogramPtr  m_Sinogram;
    //  SinogramPtr  m_BFSinogram;
    std::vector<AMatrixCol::Pointer>& m_TempCol;
    RealVolumeType::Pointer m_ErrorSino;
    std::vector<AMatrixCol::Pointer>& m_VoxelLineResponse;//CHANGED! Put a & here
    BFForwardModel* m_ForwardModel;
    UInt8Image_t::Pointer m_Mask;
    RealImageType::Pointer m_MagUpdateMap;//Hold the magnitude of the reconstuction along each voxel line
    UInt8Image_t::Pointer m_MagUpdateMask;
    //unsigned int m_VoxelUpdateType;

    int m_ZeroCount;
    Real_t m_CurrentVoxelValue;
#if ROI
    //variables used to stop the process
    Real_t* m_AverageUpdate;
    Real_t* m_AverageMagnitudeOfRecon;
#endif
    unsigned int m_ZeroSkipping;

    QGGMRF::QGGMRF_Values* m_QggmrfValues;

    //if 1 then this is NOT outside the support region; If 0 then that pixel should not be considered
    uint8_t m_BoundaryFlag[27];

    //Real_t m_HammingWindow[5][5];
    Real_t m_Theta1;
    Real_t m_Theta2;
    Real_t m_Filter[27];
    Real_t m_Neighborhood[27];

    //struct List* m_VoxelUpdateList;
    VoxelUpdateList::Pointer m_VoxelUpdateList;
};

#endif /* UPDATEYSLICE_H_ */
