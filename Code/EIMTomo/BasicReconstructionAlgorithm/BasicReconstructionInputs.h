#ifndef COMPUTATIONINPUTS_H_
#define COMPUTATIONINPUTS_H_

#include <stdio.h> //For all other declarations int,FILE etc
#include <string>

#include "BasicReconstructionStructures.h"

class BasicReconstructionInputs
{
  public:
    BasicReconstructionInputs();
    virtual ~BasicReconstructionInputs();

    int CI_ParseInput(int argc, char **argv, CommandLineInputs* Input);
    void CI_ReadParameterFile(FILE *Fp, CommandLineInputs* ParsedInput, Sino* Sinogram, Geom* Geometry);
    void CI_InitializeSinoParameters(Sino* Sinogram, CommandLineInputs* ParsedInput);
    void CI_InitializeGeomParameters(Sino* Sinogram, Geom* Geometry, CommandLineInputs* ParsedInput);

    /**
     * @brief Copys the std::string contents into a newly malloc'ed char array which
     * the programmer will need to free when they are finished with it.
     * @param fname The filename to copy
     */
    char* copyFilenameToNewCharBuffer(const std::string &fname);

  private:

    BasicReconstructionInputs(const BasicReconstructionInputs&); // Copy Constructor Not Implemented
    void operator=(const BasicReconstructionInputs&); // Operator '=' Not Implemented

};

#endif //ComputationInputs
