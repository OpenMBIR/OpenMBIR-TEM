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

	//Compute A^t * A
	
	//Compute T1=(A^t*A)^-1 
	
	//Compute T2=A^t*Y
	
	//Compute T1*T2 (2 X 2 matrix with 2 X 1 vector)

  setErrorCondition(0);
  setErrorMessage("");
  notify("Done ComputeGainsOffsets", 0, UpdateProgressMessage);
}
