/*
 * TomoFilter.cpp
 *
 *  Created on: Dec 6, 2011
 *      Author: mjackson
 */

#include "TomoFilter.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoFilter::TomoFilter() :
m_Inputs(NULL),
m_Sinogram(NULL),
m_Geometry(NULL)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoFilter::~TomoFilter()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoFilter::execute()
{
  setErrorCondition(-1);
  setErrorMessage("TomoFilter must be subclassed to be usable");
  notify(getErrorMessage().c_str(), 100, UpdateProgressMessage);
}
