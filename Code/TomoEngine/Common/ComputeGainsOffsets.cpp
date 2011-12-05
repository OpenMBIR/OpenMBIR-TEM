/*
 * ComputeGainsOffets.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: mjackson
 */

#include "ComputeGainsOffsets.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeGainsOffsets::ComputeGainsOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeGainsOffsets::~ComputeGainsOffsets()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeGainsOffsets::execute()
{
  setErrorCondition(0);
  setErrorMessage("");
  notify("Done ComputeGainsOffsets", 0, UpdateProgressMessage);


}
