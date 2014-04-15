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

#include "FilterPipeline.h"

#include "MXA/MXA.h"

#if MXA_BENCHMARKS
#define START_CLOCK()\
  unsigned long long int millis;\
  millis = MXA::getMilliSeconds();
#else
#define START_CLOCK() unsigned long long int millis = 0;\
  millis = 0;
#endif


#define CHECK_FOR_CANCELED(FuncClass, Message, name)\
    if (this->getCancel() ) { \
              pipelineProgressMessage(#Message);\
              pipelineProgress(0);\
              pipelineFinished();\
              return;}\


#define CHECK_FOR_ERROR(FuncClass, Message, err)\
    if(err < 0) {\
      setErrorCondition(err);\
      std::string msg = std::string(Message);\
      pipelineErrorMessage(msg.c_str());\
      pipelineProgress(0);\
      pipelineFinished();\
      return;   }



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FilterPipeline::FilterPipeline() :
    Observer(),
    m_ErrorCondition(0),
    m_Cancel(false)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FilterPipeline::~FilterPipeline()
{
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FilterPipeline::pipelineFinished()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FilterPipeline::setCancel(bool value)
{
  this->m_Cancel = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool FilterPipeline::getCancel()
{
  return m_Cancel;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FilterPipeline::run()
{
  execute();
  pipelineFinished();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FilterPipeline::execute()
{
    //int err = 0;

    std::stringstream ss;


    ss.str("");
    ss << "FilterPipeline Complete";
    pipelineProgressMessage(ss.str());
}
