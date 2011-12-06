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
#ifndef _AbstractFilter_H_
#define _AbstractFilter_H_


#include "MXA/MXA.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "TomoEngine/Common/Observable.h"


/**
 * @class AbstractFilter AbstractFilter.h MXA/Common/AbstractFilter.h
 * @brief This class is the basic class to subclass when creating a new Filter for
 * DREAM.3D. The subclass must implement at least the  execute method. If an
 * error occurs during the execution of the filter set the errorCondition to
 * a non zero value and optionally use the setErrorMessage() method to explain what the
 * error was. This class also inherits from Observable so that the filter can
 * provide updates to the user interface during execution.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Nov 28, 2011
 * @version 1.0
 */
class MXA_EXPORT AbstractFilter : public Observable
{
  public:
    MXA_SHARED_POINTERS(AbstractFilter)
    MXA_STATIC_NEW_MACRO(AbstractFilter)
    MXA_TYPE_MACRO_SUPER(AbstractFilter, Observable)

    virtual ~AbstractFilter();


    MXA_INSTANCE_PROPERTY(int, ErrorCondition);

    MXA_INSTANCE_STRING_PROPERTY(ErrorMessage);

    /**
     * @brief This method should be fully implemented in subclasses.
     */
    virtual void execute();

  protected:
    AbstractFilter();

  private:
    AbstractFilter(const AbstractFilter&); // Copy Constructor Not Implemented
    void operator=(const AbstractFilter&); // Operator '=' Not Implemented
};




#endif /* _AbstractFilter_H_  */
