#ifndef COMPUTATIONINPUTS_H_
#define COMPUTATIONINPUTS_H_

#include <stdio.h>
#include <string>
#include "ScaleOffsetStructures.h"




class ScaleOffsetCorrectionParser
{
  public:
    ScaleOffsetCorrectionParser();
    virtual ~ScaleOffsetCorrectionParser();

    int parseArguments(int argc, char **argv, TomoInputs* Input);
    void readParameterFile(FILE *Fp,TomoInputs* ParsedInput,Sino* Sinogram,Geom* Geometry);
    void initializeSinoParameters(Sino* Sinogram,TomoInputs* ParsedInput);
    void initializeGeomParameters(Sino* Sinogram,Geom* Geometry,TomoInputs* ParsedInput);
    DATA_TYPE absMaxArray(DATA_TYPE* Array, uint16_t NumElts);

    /**
     * @brief Copys the std::string contents into a newly malloc'ed char array which
     * the programmer will need to free when they are finished with it.
     * @param fname The filename to copy
     */
    char* copyFilenameToNewCharBuffer( const std::string &fname);


  private:
    uint64_t startm;
    uint64_t stopm;

    ScaleOffsetCorrectionParser(const ScaleOffsetCorrectionParser&); // Copy Constructor Not Implemented
    void operator=(const ScaleOffsetCorrectionParser&); // Operator '=' Not Implemented

};

#endif //ComputationInputs
