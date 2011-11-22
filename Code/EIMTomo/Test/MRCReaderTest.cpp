/*
 * MRCReaderTest.cpp
 *
 *  Created on: Nov 22, 2011
 *      Author: mjackson
 */

#include <stdlib.h>

#include <string>
#include <iostream>

#include "EIMTomo/IO/MRCReader.h"
#include "EIMTomo/IO/MRCHeader.h"

const std::string filepath("/Users/Shared/Data/TomographyData/TiO2Ps100kRun3/Run3TiO2PS100k_2.ali");

int main(int argc, char **argv)
{
  std::cout << "Testing file \n  " << filepath << std::endl;

  MRCReader reader(false);
  int err = reader.read(filepath);
  std::cout << "Error Code from Reading: " << err << std::endl;

  if (err < 0)
  {
    return EXIT_FAILURE;
  }

  void* data = reader.getDataPointer();

  free(data);

  std::cout << "Done Reading File" << std::endl;
  return EXIT_SUCCESS;
}


