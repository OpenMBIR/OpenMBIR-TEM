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

#include "VoxelUpdate.h"

#include <stdio.h>

#include "TomoEngine/Common/EIMTime.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCConstants.h"
#include "TomoEngine/Common/DerivOfCostFunc.hpp"

#define START startm = EIMTOMO_getMilliSeconds();
#define STOP stopm = EIMTOMO_getMilliSeconds();
#define PRINTTIME printf( "%6.3f seconds used by the processor.\n", ((double)stopm-startm)/1000.0);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdate::VoxelUpdate() :
    TomoFilter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdate::~VoxelUpdate()
{
}

/*****************************************************************************
 //Finds the min and max of the neighborhood . This is required prior to calling
 solve()
 *****************************************************************************/
void VoxelUpdate::minMax(DATA_TYPE *low, DATA_TYPE *high)
{
  DATA_TYPE THETA1 = m_VoxelUpdateValues->THETA1;
  DATA_TYPE THETA2 = m_VoxelUpdateValues->THETA2;
  DATA_TYPE V = m_VoxelUpdateValues->V;

  *low = m_VoxelUpdateValues->NEIGHBORHOOD[0][0][0];
  *high = m_VoxelUpdateValues->NEIGHBORHOOD[0][0][0];

  for (int32_t i = 0; i < 3; i++)
  {
    for (int32_t j = 0; j < 3; j++)
    {
      for (int32_t k = 0; k < 3; k++)
      {
        //  if(NEIGHBORHOOD[i][j][k] != 0)
        //  printf("%lf ", NEIGHBORHOOD[i][j][k]);

        if(m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k] < *low)
        {
          *low = m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k];
        }
        if(m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k] > *high)
        {
          *high = m_VoxelUpdateValues->NEIGHBORHOOD[i][j][k];
        }
      }
      //  printf("\n");
    }
  }

  if(THETA2 != 0)
  {
    *low = (*low > (V - (THETA1 / THETA2)) ? (V - (THETA1 / THETA2)) : *low);
    *high = (*high < (V - (THETA1 / THETA2)) ? (V - (THETA1 / THETA2)) : *high);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdate::execute()
{}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
