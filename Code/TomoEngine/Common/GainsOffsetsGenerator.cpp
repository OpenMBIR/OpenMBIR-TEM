/*
 * GainsOffsetsGenerator.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "GainsOffsetsGenerator.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsGenerator::GainsOffsetsGenerator()
{

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsGenerator::~GainsOffsetsGenerator()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsGenerator::execute()
{




  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Generating the Gains and Offsets Data", 0, UpdateProgressMessage);

}
