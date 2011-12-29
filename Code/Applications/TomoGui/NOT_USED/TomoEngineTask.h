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

#ifndef _TOMOENGINE_TASK_H_
#define _TOMOENGINE_TASK_H_


#include <QtCore/QObject>
#include <QtCore/QThread>
#include <QtCore/QString>
#include <QtGui/QImage>


#include "MXA/Common/MXASetGetMacros.h"

#include "QtSupport/ProcessQueueTask.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/SOC/SOCEngine.h"
#include "TomoEngine/SOC/SOCStructures.h"

#define UPDATE_PROGRESS(m, p)\
  emit progressTextChanged( (m) );\
  emit progressValueChanged( (p) );

/**
* @class TomoEngineTask TomoEngineTask.h EmMpm/GUI/TomoEngineTask.h
* @brief This is the wrapper code for the code. This is called as a "worker" class from a separate thread
* of execution in order to not lock up the GUI.
* @author Michael A. Jackson for BlueQuartz Software
* @date Dec 20, 2009
* @version 1.0
*/
class TomoEngineTask : public ProcessQueueTask
{

  Q_OBJECT;

  public:
    TomoEngineTask(QObject* parent = 0);
    virtual ~TomoEngineTask();

    MXA_INSTANCE_PROPERTY(TomoInputsPtr, TomoInputs)
    MXA_INSTANCE_PROPERTY(SinogramPtr, Sinogram)
    MXA_INSTANCE_PROPERTY(GeometryPtr, Geometry)
    MXA_INSTANCE_PROPERTY(ScaleOffsetParamsPtr, NuisanceParams)

    virtual void run();

  public slots:

    /**
     * @brief Slot to receive a signal to cancel the operation
     */
    void cancel();

  signals:
     /**
     * @brief Signal sent when the encoder has a message to relay to the GUI or other output device.
     */
      void progressTextChanged (QString progressText );



  private:
      SOCEngine::Pointer m_Engine;

    TomoEngineTask(const TomoEngineTask&); // Copy Constructor Not Implemented
    void operator=(const TomoEngineTask&); // Operator '=' Not Implemented

};





#endif /* _TOMOENGINE_TASK_H_ */
