/*
 * JointEstimationCLIParser.h
 *
 *  Created on: Sep 19, 2011
 *      Author: mjackson
 */

#ifndef JOINTESTIMATIONCLIPARSER_H_
#define JOINTESTIMATIONCLIPARSER_H_

#include <string>

#include "MotionCorrectionInputs.h"


/*
 *
 */
class JointEstimationCLIParser
{
  public:
    JointEstimationCLIParser();
    virtual ~JointEstimationCLIParser();

    /**
     * @brief
     * @param argc
     * @param argv
     * @param inputs
     * @return
     */
    int parseCLIArguments(int argc,char *argv[], CommandLineInputs* inputs);

    /**
     * @brief Copys the std::string contents into a newly malloc'ed char array which
     * the programmer will need to free when they are finished with it.
     * @param fname The filename to copy
     */
    char* copyFilenameToNewCharBuffer( const std::string &fname);

  private:
    JointEstimationCLIParser(const JointEstimationCLIParser&);    // Copy Constructor Not Implemented
      void operator=(const JointEstimationCLIParser&);  // Operator '=' Not Implemented
};

#endif /* JOINTESTIMATIONCLIPARSER_H_ */
