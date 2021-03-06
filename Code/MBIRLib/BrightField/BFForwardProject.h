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



#ifndef FORWARDPROJECT_H_
#define FORWARDPROJECT_H_

#include <iostream>


#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/Common/Observable.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"

#include "MBIRLib/Common/AMatrixCol.h"
#include "BFForwardModel.h"


/**
 * @class BFForwardProject BFForwardProject.h TomoEngine/SOC/BFForwardProject.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @author Singanallur Venkatakrishnan (Purdue University)
 * @date Dec 12, 2011
 * @version 1.0
 */
class MBIRLib_EXPORT BFForwardProject
{
  public:
    BFForwardProject(Sinogram* sinogram,
                     Geometry* geometry,
                     std::vector<AMatrixCol::Pointer>& tempCol,
                     std::vector<AMatrixCol::Pointer>& voxelLineResponse,
                     RealVolumeType::Pointer yEst,
                     BFForwardModel* forwardModel,
                     uint16_t tilt,
                     Observable* obs);

    virtual ~BFForwardProject();

    void operator()() const;

  private:
    Sinogram* m_Sinogram;
    Geometry* m_Geometry;
    std::vector<AMatrixCol::Pointer> m_TempCol;
    std::vector<AMatrixCol::Pointer> m_VoxelLineResponse;
    RealVolumeType::Pointer m_YEstimate;
    BFForwardModel* m_ForwardModel;
    uint16_t m_Tilt;
    Observable* m_Observable;
};

#endif /* FORWARDPROJECT_H_ */
