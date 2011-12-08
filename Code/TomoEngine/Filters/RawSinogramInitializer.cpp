/*
 * RawSinogramInitializer.cpp
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#include "RawSinogramInitializer.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RawSinogramInitializer::RawSinogramInitializer()
{

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RawSinogramInitializer::~RawSinogramInitializer()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RawSinogramInitializer::execute()
{
  // If an error occurs, clean up any memory, call "setErrorCondition(-1)" and
  // also setErrorMessage("Something went wrong"); and then return





  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the Raw Input file", 0, UpdateProgressMessage);

}
