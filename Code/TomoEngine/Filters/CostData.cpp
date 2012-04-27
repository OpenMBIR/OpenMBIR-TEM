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

#include "CostData.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CostData::CostData() :
m_File(NULL)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CostData::~CostData()
{
  if (NULL != m_File)
  {
    fclose(m_File);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int CostData::numberOfCosts()
{
  return m_Cost.size();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CostData::printCosts(std::ostream &out)
{
  out << "Cost Values -----------------------------" << std::endl;
  for(std::vector<Real_t>::size_type ii = 0 ; ii < m_Cost.size(); ii++)
  {
    out << m_Cost[ii] << std::endl;
  }
  out << "-----------------------------------------" << std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int CostData::initOutputFile(const std::string &filepath)
{
  m_File = fopen(filepath.c_str(), "wb");
  if (m_File != NULL)
    return 1;
  else
    return -1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int CostData::addCostValue(Real_t value)
{
  m_Cost.push_back(value);
  if(m_Cost[m_Cost.size() - 1] - m_Cost[m_Cost.size() - 2] > 0)
  {
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int CostData::writeCostValue(Real_t value)
{
  if (m_File != NULL) {
    fwrite(&(value), sizeof(Real_t), 1, m_File);
    return 1;
  }
  return -1;
}
