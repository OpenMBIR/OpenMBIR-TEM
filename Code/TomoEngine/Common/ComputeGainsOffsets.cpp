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
  // If an error occurs, clean up any memory, call "setErrorCondition(-1)" and
  // also setErrorMessage("Something went wrong"); and then return


  setErrorCondition(0);
  setErrorMessage("");
  notify("Done ComputeGainsOffsets", 0, UpdateProgressMessage);
}
