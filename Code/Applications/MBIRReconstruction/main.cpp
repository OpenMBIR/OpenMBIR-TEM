/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "TomoEngine/TomoEngine.h"

#include <stdlib.h>
#include <string.h>

#include <string>
#include <iostream>


// MXA Includes
#include "MXA/Utilities/MXADir.h"
#include "MXA/Utilities/MXAFileInfo.h"

#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/SOC/SOCStructures.h"
#include "TomoEngine/SOC/SOCEngine.h"
#include "TomoEngine/SOC/MultiResolutionSOC.h"
#include "MBIRReconstructionArgsParser.h"

int main(int argc, char **argv)
{
    std::cout << "Starting MBIR Reconstruction Version " << TomoEngine::Version::Complete << std::endl;



    MultiResolutionSOC::Pointer engine = MultiResolutionSOC::New();

    MBIRReconstructionArgsParser argParser;
    int err = argParser.parseArguments(argc, argv, engine);
    if(err < 0)
    {
        std::cout << "Error Parsing the arguments." << std::endl;
        return EXIT_FAILURE;
    }


#if 0
    char path1[MAXPATHLEN]; // This is a buffer for the text
    ::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
    getcwd(path1, MAXPATHLEN);
    std::cout << "Current Working Directory: " << path1 << std::endl;
#endif
    // Make sure the output directory is created if it does not exist
    // Get the absolute path to the users choice for the output file
    std::string parentPath = MXAFileInfo::absolutePath(engine->getOutputFile());


    if(MXADir::exists(parentPath) == false)
    {
        std::cout << "Output Directory '" << parentPath << "' does NOT exist. Attempting to create it." << std::endl;
        if(MXADir::mkdir( parentPath, true) == false)
        {
            std::cout << "Error creating the output directory '" << parentPath << "'\n   Exiting Now." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "Output Directory Created." << std::endl;
    }
#if 0
    ::memset(path1, 0, MAXPATHLEN); // Initialize the string to all zeros.
    getcwd(path1, MAXPATHLEN);
    std::cout << "Current Working Directory: " << path1 << std::endl;
#endif

    std::stringstream ss;


    // Run the reconstruction
    engine->execute();
    if(engine->getErrorCondition() < 0)
    {
        std::cout << "Error Reconstructing the Data" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Completed MBIR Reconstruction Run" << std::endl;

    return EXIT_SUCCESS;
}

