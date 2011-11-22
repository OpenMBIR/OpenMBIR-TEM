
#include "EIMTomo/EIMTomoConfiguration.h"


// C Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#ifdef CMP_HAVE_SYS_PARAM_H
#include <sys/param.h>  // MAXPATHLEN definition
#else
#error sys/paramh is needed for this code and was not found on your system
#endif
// C++ includes
#include <string>
#include <iostream>

// EIMTomo Includes
#include "EIMTomo/EIMTomo.h"
#include "EIMTomo/common/EIMTime.h"
#include "EIMTomo/common/allocate.h"

#include "EIMTomo/IO/VTKFileWriters.hpp"
#include "EIMTomo/IO/RawGeometryWriter.h"

// MXA Includes
#include "MXA/Utilities/MXADir.h"
#include "ScaleOffsetMotionCorrectionConstants.h"
#include "ScaleOffsetMotionCorrectionParser.h"
#include "ScaleOffsetMotionCorrectionEngine.h"

#define START startm = EIMTOMO_getMilliSeconds();
#define STOP stopm = EIMTOMO_getMilliSeconds();
#define PRINTTIME printf( "%6.3f seconds used by the processor.\n", ((double)stopm-startm)/1000.0);

/*double
 CE => Computation Engine
 CI => Computation Inputs
 N_ => Number of ..
 
 */
#define MAKE_OUTPUT_FILE(Fp, err, outdir, filename)\
{\
std::string filepath(outdir);\
filepath =+ MXADir::Separator;\
filepath.append(filename);\
Fp = fopen(filepath.c_str(),"wb");\
if (Fp == NULL) { std::cout << "Error Opening Output file " << filepath << std::endl; err = 1; }\
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
	int32_t error;
	//FILE* Fp = NULL;
	TomoInputs inputs;
	Sino sinogram;
	Geom geometry;
	uint64_t startm;
	uint64_t stopm;
	
	START;
	
	error = -1;
	ScaleOffsetCorrectionParser soci;
	
	error = soci.parseArguments(argc, argv, &inputs);
	if(error != -1)
	{
		printf("***************\nParameter File: %s\nSinogram File: %s\nInitial Reconstruction File: %s\nOutput Directory: %s\nOutput File: %s\n***************\n",
			   inputs.ParamFile.c_str(),
			   inputs.SinoFile.c_str(),
			   inputs.InitialRecon.c_str(),
			   inputs.outputDir.c_str(),
			   inputs.OutputFile.c_str());
		
		char path1[MAXPATHLEN];  // This is a buffer for the text
		::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
		getcwd(path1, MAXPATHLEN);
		std::cout << "Current Working Directory: " << path1 << std::endl;
		
		
		// Make sure the output directory is created if it does not exist
		if(MXADir::exists(inputs.outputDir) == false)
		{
			std::cout << "Output Directory '" << inputs.outputDir << "' does NOT exist. Attempting to create it." << std::endl;
			if(MXADir::mkdir(inputs.outputDir, true) == false)
			{
				std::cout << "Error creating the output directory '" << inputs.outputDir << "'\n   Exiting Now." << std::endl;
				return EXIT_FAILURE;
			}
			std::cout << "Output Directory Created." << std::endl;
		}
		::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
		getcwd(path1, MAXPATHLEN);
		std::cout << "Current Working Directory: " << path1 << std::endl;
		
		error = soci.readParameterFile(inputs.ParamFile, &inputs, &sinogram, &geometry);
		if (error < 0)
		{
			std::cout << "Error Opening Parameter file and parsing the data within." << std::endl;
			return EXIT_FAILURE;
		}
		//Based on the inputs , calculate the "other" variables in the structure definition
		soci.initializeSinoParameters(&sinogram, &inputs);
		//CI_MaskSinogram(&OriginalSinogram,&MaskedSinogram);
		soci.initializeGeomParameters(&sinogram, &geometry, &inputs);
	}
	
	// Run the reconstruction
	SOCEngine soce(&sinogram, &geometry, &inputs);
	error = soce.mapicdReconstruct();
	if(error < 0)
	{
		std::cout << "Error (" << error << ") During Reconstruction" << std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "Final Dimensions of Object: " << std::endl;
	std::cout << "  Nx = " << geometry.N_x << std::endl;
	std::cout << "  Ny = " << geometry.N_y << std::endl;
	std::cout << "  Nz = " << geometry.N_z << std::endl;
	
	// Write the vtk output file
	VTKRectilinearGridFileWriter vtkWriter;
	vtkWriter.setWriteBinaryFiles(true);
	DimsAndRes dimsAndRes;
	dimsAndRes.xpoints = geometry.N_x;
	dimsAndRes.ypoints = geometry.N_y;
	dimsAndRes.zpoints = geometry.N_z;
	dimsAndRes.resx = 1.0f;
	dimsAndRes.resy = 1.0f;
	dimsAndRes.resz = 1.0f;
	std::string vtkFile(inputs.outputDir);
	vtkFile = vtkFile.append(MXADir::getSeparator()).append("Output.vtk");
	
	VtkScalarWriter* w0 = static_cast<VtkScalarWriter*>(new TomoOutputScalarWriter(&geometry));
	std::vector<VtkScalarWriter*> scalarsToWrite;
	w0->m_WriteBinaryFiles = true;
	scalarsToWrite.push_back(w0);
	
	error = vtkWriter.write<DimsAndRes>(vtkFile, &dimsAndRes, scalarsToWrite);
	if (error < 0)
	{
		std::cout << "Error writing vtk file '" << vtkFile << "'" << std::endl;
	}
	
	// Write a raw binary output file
	RawGeometryWriter binWriter(&geometry);
	error = binWriter.writeFile(inputs.OutputFile);
	if (error < 0)
	{
		std::cout << "Error writing Raw binary file '" << inputs.OutputFile << "'" << std::endl;
	}
	
	
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

