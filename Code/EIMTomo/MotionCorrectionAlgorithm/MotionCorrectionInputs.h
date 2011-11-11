#ifndef MOTION_CORRECTION_COMPUTATIONINPUTS_H_
#define MOTION_CORRECTION_COMPUTATIONINPUTS_H_

#include <stdio.h> //For all other declarations int,FILE etc
#include <string>
#include "MotionCorrectionStructures.h"

#define START_SLICE 0 //This is used to select whcih slice to reconstruct - REMOVE this in the final versions of the code
#define END_SLICE 3//This is used to select whcih slice to reconstruct - REMOVE this in the final versions of the code



class MotionCorrectionInputs
{

  public:
    MotionCorrectionInputs();
    virtual ~MotionCorrectionInputs();

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
    MotionCorrectionInputs(const MotionCorrectionInputs&); // Copy Constructor Not Implemented
    void operator=(const MotionCorrectionInputs&); // Operator '=' Not Implemented
};

#endif //ComputationInputs
