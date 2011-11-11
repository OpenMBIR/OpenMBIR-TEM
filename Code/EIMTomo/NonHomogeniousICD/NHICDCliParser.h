/*
 * NHICDCliParser.h
 *
 *  Created on: Sep 19, 2011
 *      Author: mjackson
 */

#ifndef NHICDCLIPARSER_H_
#define NHICDCLIPARSER_H_

#include <string>

#include "NHICDInputs.h"

/*
 *
 */
class NHICDCliParser
{
  public:
    NHICDCliParser();
    virtual ~NHICDCliParser();

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
   NHICDCliParser(const NHICDCliParser&);    // Copy Constructor Not Implemented
    void operator=(const NHICDCliParser&);  // Operator '=' Not Implemented
};

#endif /* NHICDCLIPARSER_H_ */
