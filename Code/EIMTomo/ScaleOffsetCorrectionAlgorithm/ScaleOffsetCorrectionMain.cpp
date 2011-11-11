#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "EIMTomo/EIMTomo.h"
#include "EIMTomo/common/EIMTime.h"
#include "EIMTomo/common/allocate.h"

#include "ScaleOffsetCorrectionInputs.h"
#include "ScaleOffsetCorrectionEngine.h"

#define START startm = EIMTOMO_getMilliSeconds();
#define STOP stopm = EIMTOMO_getMilliSeconds();
#define PRINTTIME printf( "%6.3f seconds used by the processor.\n", ((double)stopm-startm)/1000.0);

/*double
 CE => Computation Engine
 CI => Computation Inputs
 N_ => Number of ..

 */

int main(int argc, char** argv)
{
  int16_t error, i, j, k;
  FILE* Fp;
  CommandLineInputs ParsedInput;
  Sino Sinogram;
  Geom Geometry;
  DATA_TYPE *buffer = (DATA_TYPE*)get_spc(1, sizeof(DATA_TYPE));
  uint64_t startm;
  uint64_t stopm;

  START;

  error = -1;
  ScaleOffsetCorrectionInputs soci;

  error = soci.CI_ParseInput(argc, argv, &ParsedInput);
  if(error != -1)
  {
    printf("%s %s %s %s\n", ParsedInput.ParamFile, ParsedInput.SinoFile, ParsedInput.InitialRecon, ParsedInput.OutputFile);
    //Read the paramters into the structures
    Fp = fopen(ParsedInput.ParamFile, "r");
    soci.CI_ReadParameterFile(Fp, &ParsedInput, &Sinogram, &Geometry);
    fclose(Fp);
    //Based on the inputs , calculate the "other" variables in the structure definition
    soci.CI_InitializeSinoParameters(&Sinogram, &ParsedInput);
    //CI_MaskSinogram(&OriginalSinogram,&MaskedSinogram);
    soci.CI_InitializeGeomParameters(&Sinogram, &Geometry, &ParsedInput);
  }

  SOCEngine soce;
  error = soce.CE_MAPICDReconstruct(&Sinogram, &Geometry, &ParsedInput);
  Fp = fopen(ParsedInput.OutputFile, "w");

  printf("Main\n");
  printf("Final Dimensions of Object Nz=%d Nx=%d Ny=%d\n", Geometry.N_z, Geometry.N_x, Geometry.N_y);

  for (i = 0; i < Geometry.N_y; i++)
  {
    for (j = 0; j < Geometry.N_x; j++)
      for (k = 0; k < Geometry.N_z; k++)
      {
        buffer = &Geometry.Object[k][j][i];
        fwrite(buffer, sizeof(DATA_TYPE), 1, Fp);
      }
    printf("%d\n", i);
  }

  fclose(Fp);
  STOP;
  PRINTTIME;
  // if(error < 0)
  // {
  // //TODO:clean up memory here
  // return EXIT_FAILURE;
  // }
  // //TODO:free memory
  //
  // return EXIT_SUCCESS;
  return 0;
}

