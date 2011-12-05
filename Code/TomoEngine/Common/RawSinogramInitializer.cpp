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
RawSinogramInitializer::RawSinogramInitializer() :
m_Inputs(NULL),
m_Sinogram(NULL)
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




  setErrorCondition(0);
  setErrorMessage("");
  notify("Done Reading the Raw Input file", 0, UpdateProgressMessage);

}
