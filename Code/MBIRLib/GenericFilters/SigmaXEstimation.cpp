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



#include "SigmaXEstimation.h"

#include <limits>

#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/IOFilters/MRCHeader.h"
#include "MBIRLib/IOFilters/MRCReader.h"
#include "MBIRLib/Reconstruction/ReconstructionConstants.h"


namespace Detail
{
    const float DegToRad = 0.017453292519943f;

template<typename T>
void calcMinMax(T* data, int total, Real_t &min, Real_t &max, Real_t &sum2)
{
  for (int i = 0; i < total; i++)
  {
      if(data[i] > max) max = data[i];
      if(data[i] < min) min = data[i];
      sum2 += (data[i]);
  }
}

template<typename T>
void calcAvgDeviation(T* data, int total, Real_t &dev)
{
    dev=0;
    Real_t mean=0;
    for (int i = 0; i < total; i++)
    {
        mean+=data[i];
    }
    mean/=total;
    for (int i = 0; i < total; i++)
    {
        dev+=fabs(data[i]-mean);
    }
    dev/=total;
}

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SigmaXEstimation::SigmaXEstimation() :
m_SampleThickness(100.0),
m_TiltAngles(0),
m_DefaultOffset(0.0),
m_TargetGain(1.0),
m_SigmaXEstimate(0.0)
{
    m_XDims[0] = -1;
    m_XDims[1] = -1;
    m_YDims[0] = -1;
    m_YDims[1] = -1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SigmaXEstimation::~SigmaXEstimation()
{

}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SigmaXEstimation::execute()
{

    MRCHeader header;
    ::memset(&header, 0, 1024);
    header.feiHeaders = NULL;
    MRCReader::Pointer reader = MRCReader::New(true);
    int err = reader->readHeader(m_InputFile, &header);
    if (err < 0)
    {
        FREE_FEI_HEADERS( header.feiHeaders )
        notify("Error reading the MRC input file", 0, Observable::UpdateErrorMessage);
        setErrorCondition(-1);
        return;
    }

  //Make sure min/max have been set otherwise just use the dimensions of the data from the header
    if (m_XDims[0] < 0) { m_XDims[0] = 0; }
    if (m_XDims[1] < 0) { m_XDims[1] = header.nx;}
    if (m_YDims[0] < 0) { m_YDims[0] = 0;}
    if (m_YDims[1] < 0) { m_YDims[1] = header.ny;}

    int voxelMin[3] = { m_XDims[0], m_YDims[0], 0};
    int voxelMax[3] = { m_XDims[1], m_YDims[1], 0};
    Real_t sum1 = 0;
//    Real_t targetMin = std::numeric_limits<Real_t>::max();
 //   Real_t targetMax = std::numeric_limits<Real_t>::min();
//    Real_t min = std::numeric_limits<Real_t>::max();
//    Real_t max = std::numeric_limits<Real_t>::min();

    std::vector<Real_t> sum2s(header.nz);
    int voxelCount = (m_XDims[1] - m_XDims[0]) * (m_YDims[1] - m_YDims[0]);
    float progress = 0.0;
  // Loop over each tilt angle to compute the Target Gain Estimation
    for(int i_theta = 0; i_theta < header.nz; ++i_theta)
    {
        voxelMin[2] = i_theta;
        voxelMax[2] = i_theta;
        err = reader->read(m_InputFile, voxelMin, voxelMax); // This should read the subvolume define by voxelMin and voxelMax
        if(err < 0)
        {
            break;
        }
        progress = (i_theta/header.nz) * 100.0f;


        Real_t sum2 = 0;
        switch(header.mode)
        {
        case 0:
            //calcMinMax<uint8_t>(static_cast<uint8_t*>(reader->getDataPointer()), voxelCount, min, max, sum2);
            Detail::calcAvgDeviation<uint8_t>(static_cast<uint8_t*>(reader->getDataPointer()), voxelCount, sum2);
            break;
        case 1:
            //calcMinMax<int16_t>(static_cast<int16_t*>(reader->getDataPointer()), voxelCount, min, max, sum2);
            Detail::calcAvgDeviation<int16_t>(static_cast<int16_t*>(reader->getDataPointer()), voxelCount, sum2);
            break;
        case 2:
            //calcMinMax<float>(static_cast<float*>(reader->getDataPointer()), voxelCount, min, max, sum2);
            Detail::calcAvgDeviation<float>(static_cast<float*>(reader->getDataPointer()), voxelCount, sum2);
            break;
        case 3:
            break;
        case 4:
            break;
        case 6:
            //calcMinMax<uint16_t>(static_cast<uint16_t*>(reader->getDataPointer()), voxelCount, min, max, sum2);
            Detail::calcAvgDeviation<uint16_t>(static_cast<uint16_t*>(reader->getDataPointer()), voxelCount, sum2);
            break;
        case 16:
            break;
        }
//        if (min < targetMin) { targetMin = min; }
//        if (max > targetMax) { targetMax = max; }
        sum2s[i_theta] = sum2;

        notify("Estimating Target Gain and Sigma X from Data. ", (int)progress, Observable::UpdateProgressValueAndMessage);
    }


    //modify it based on any knowledge they have about the tx. attenuation

    // Now Calculate the Sigma X estimation
    Real_t cosine = 0.0;
    for(int i_theta = 0; i_theta < header.nz; ++i_theta)
    {
        //Subtract off any offset in the data
  //      sum2s[i_theta] -= min * voxelCount;

  //      sum2s[i_theta] /= voxelCount * m_TargetGain;
        sum2s[i_theta]/=m_TargetGain;

        if(m_TiltAngles == 0)
        {
            cosine = cos(header.feiHeaders[i_theta].a_tilt * (Detail::DegToRad));
        }
        else
        {
            cosine = cos(header.feiHeaders[i_theta].b_tilt * (Detail::DegToRad));
        }

        sum1 += (sum2s[i_theta] * cosine) / (m_SampleThickness);
    }

    m_SigmaXEstimate = sum1/header.nz;///10.0;
    FREE_FEI_HEADERS( header.feiHeaders );

    //  std::cout << "Estimated Target Gain: " << m_TargetGainEstimate << std::endl;
    //  std::cout << "Estimated Sigma X: " << m_SigmaXEstimate << std::endl;
    notify("Estimating Target Gain and Sigma X Complete ", 100, Observable::UpdateProgressValueAndMessage);

}
