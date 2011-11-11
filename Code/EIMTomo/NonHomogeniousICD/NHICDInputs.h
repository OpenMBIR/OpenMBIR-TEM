#ifndef NCICD_COMPUTATIONINPUTS_H_
#define NCICD_COMPUTATIONINPUTS_H_

#include <stdio.h>
#include <string>
#include "NHICDStructures.h"



class NHICDInputs
{
  public:
    NHICDInputs();
    virtual ~NHICDInputs();


    int CI_ParseInput(int argc, char **argv, CommandLineInputs* Input);
	void CI_ReadParameterFile(FILE *Fp,CommandLineInputs* ParsedInput,Sino* Sinogram,Geom* Geometry);
	void CI_InitializeSinoParameters(Sino* Sinogram,CommandLineInputs* ParsedInput);
	void CI_InitializeGeomParameters(Sino* Sinogram,Geom* Geometry,CommandLineInputs* ParsedInput);


  /**
   * @brief Copys the std::string contents into a newly malloc'ed char array which
   * the programmer will need to free when they are finished with it.
   * @param fname The filename to copy
   */
  char* copyFilenameToNewCharBuffer( const std::string &fname);



  private:
	NHICDInputs(const NHICDInputs&); // Copy Constructor Not Implemented
  void operator=(const NHICDInputs&); // Operator '=' Not Implemented

	};

#endif //ComputationInputs
