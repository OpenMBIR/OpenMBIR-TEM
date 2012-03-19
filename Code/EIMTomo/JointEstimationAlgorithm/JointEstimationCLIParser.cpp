/*
 * JointEstimationCLIParser.cpp
 *
 *  Created on: Sep 19, 2011
 *      Author: mjackson
 */

#include "JointEstimationCLIParser.h"

#include <string.h>

#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>
#include "TomoEngine/TomoEngineVersion.h"

JointEstimationCLIParser::JointEstimationCLIParser()
{

}

JointEstimationCLIParser::~JointEstimationCLIParser()
{

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
char* JointEstimationCLIParser::copyFilenameToNewCharBuffer(const std::string &fname)
{
  std::string::size_type size = fname.size() + 1;
  char* buf = NULL;
  if (size > 1)
  {
    buf = (char*)malloc(size);
    ::memset(buf, 0, size);
    strncpy(buf, fname.c_str(), size - 1);
  }
  return buf;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int JointEstimationCLIParser::parseCLIArguments(int argc,char *argv[], TomoInputs* Input)
{
  if ( NULL == Input)
  {
    printf("The EIMTomo_Inputs pointer was null. Returning early.\n");
    return -1;
  }

  TCLAP::CmdLine cmd("", ' ', TomoEngine::Version::Complete);

  TCLAP::ValueArg<std::string> in_inputFile("i", "inputfile", "Input Data File", true, "", "");
  cmd.add(in_inputFile);
  TCLAP::ValueArg<std::string> in_outputFile("o", "outputfile", "The Output File", true, "", "");
  cmd.add(in_outputFile);
  TCLAP::ValueArg<std::string> in_paramFile("p", "paramfile", "The Parameter File", true, "", "");
  cmd.add(in_paramFile);
  TCLAP::ValueArg<std::string> in_sinoFile("s", "sinofile", "The Sinogram File", true, "", "");
  cmd.add(in_sinoFile);
  TCLAP::ValueArg<int> in_numIter("n", "numIter", "Number of Iterations", true, 0, "0");
  cmd.add(in_numIter);
  TCLAP::ValueArg<double> in_sigmaX("l", "sigmax", "Sigma X Value", true, 1.0, "1.0");
  cmd.add(in_sigmaX);
  TCLAP::ValueArg<double> in_markov("m", "", "Markov Random Field Parameter", true, 0.0, "0.0");
  cmd.add(in_markov);

  if (argc < 2)
  {
    std::cout << "Joint Estimation Command Line Version " << cmd.getVersion() << std::endl;
    std::vector<std::string> args;
    args.push_back(argv[0]);
    args.push_back("-h");
    cmd.parse(args);
    return -1;
  }


  try
  {
   // int error = 0;
    cmd.parse(argc, argv);
    Input->InitialRecon = copyFilenameToNewCharBuffer(in_inputFile.getValue());
    Input->OutputFile = copyFilenameToNewCharBuffer(in_outputFile.getValue());
    Input->ParamFile = copyFilenameToNewCharBuffer(in_paramFile.getValue());
    Input->SinoFile = copyFilenameToNewCharBuffer(in_sinoFile.getValue());
    Input->NumIter = in_numIter.getValue();
    Input->SigmaX = in_sigmaX.getValue();
    Input->p = in_markov.getValue();
  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << " error: " << e.error() << " for arg " << e.argId() << std::endl;
    std::cout << "** Unknown Arguments. Displaying help listing instead. **" << std::endl;
    return -1;
  }
  return 0;
}

