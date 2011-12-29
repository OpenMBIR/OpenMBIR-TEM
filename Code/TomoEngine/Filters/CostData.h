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
#ifndef COSTDATA_H_
#define COSTDATA_H_


#include <stdio.h>

#include <vector>

#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/SOC/SOCStructures.h"


/*
 *
 */
class TomoEngine_EXPORT CostData
{
  public:
    MXA_SHARED_POINTERS(CostData);
    MXA_STATIC_NEW_MACRO(CostData);
    MXA_TYPE_MACRO(CostData);

    virtual ~CostData();

    int initOutputFile(const std::string &filepath);

    int writeCostValue(DATA_TYPE value);

    int addCostValue(DATA_TYPE value);

    int numberOfCosts();

    void printCosts(std::ostream &out);

  protected:
    CostData();

  private:
    std::vector<DATA_TYPE>  m_Cost;
    FILE* m_File;

    CostData(const CostData&); // Copy Constructor Not Implemented
    void operator=(const CostData&); // Operator '=' Not Implemented
};

#endif /* COSTDATA_H_ */