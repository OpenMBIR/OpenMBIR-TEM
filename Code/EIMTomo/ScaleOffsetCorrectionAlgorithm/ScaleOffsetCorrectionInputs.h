#ifndef COMPUTATIONINPUTS_H_
#define COMPUTATIONINPUTS_H_

#include <stdio.h>
#include <string>
#include "ScaleOffsetStructures.h"




class ScaleOffsetCorrectionInputs
{
  public:
    ScaleOffsetCorrectionInputs();
    virtual ~ScaleOffsetCorrectionInputs();

    int CI_ParseInput(int argc, char **argv, CommandLineInputs* Input);
    void CI_ReadParameterFile(FILE *Fp,CommandLineInputs* ParsedInput,Sino* Sinogram,Geom* Geometry);
    void CI_InitializeSinoParameters(Sino* Sinogram,CommandLineInputs* ParsedInput);
    void CI_InitializeGeomParameters(Sino* Sinogram,Geom* Geometry,CommandLineInputs* ParsedInput);
    DATA_TYPE AbsMaxArray(DATA_TYPE* Array, uint16_t NumElts);

    /**
     * @brief Copys the std::string contents into a newly malloc'ed char array which
     * the programmer will need to free when they are finished with it.
     * @param fname The filename to copy
     */
    char* copyFilenameToNewCharBuffer( const std::string &fname);


  private:
    uint64_t startm;
    uint64_t stopm;

    ScaleOffsetCorrectionInputs(const ScaleOffsetCorrectionInputs&); // Copy Constructor Not Implemented
    void operator=(const ScaleOffsetCorrectionInputs&); // Operator '=' Not Implemented

};

#endif //ComputationInputs
