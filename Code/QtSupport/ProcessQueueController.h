/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
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

#ifndef PROCESSQUEUECONTROLLER_H_
#define PROCESSQUEUECONTROLLER_H_

#include <QtCore/QThread>
#include <QtCore/QMutex>
#include <QtCore/QVector>

class ProcessQueue;

/**
 * @class ProcessQueueController ProcessQueueController.h QtSupport/ProcessQueueController.h
 * @brief A QThread derived class that uses Signals and Slots to controll a Task Queue
 * system. The maximum number of threads is based on the maximum number of virtual
 * cores/CPUs that are found using the function QThread::idealThreadCount();.
 *  This class works in conbination with the ProcessQueueTask
 *  class in which you define the code that you would like to be run on a separate
 *  thread. This design implementation ensures that Qt's signals and slots between
 *  each task are delivered and executed properly across threads. This is accmplished
 *  by the use of signals and slots to start/end the ProcessTask.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date May 26, 2010
 * @version 1.0
 */
class ProcessQueueController : public QThread
{
    Q_OBJECT;

  public:
    /**
     * @brief Standard Qt constructor.
     * @param parent The QObject parent of this object.
     * @return Properly constructed EmMpmThread Object
     */
    ProcessQueueController(QObject* parent = 0);

    /**
     * Destructor
     */
    virtual ~ProcessQueueController();

    /**
     * @brief Sets the flag to automatically delete the task when complete
     * @param deleteQueue Delete the Encoder task when encoding is complete.
     */
    void setAutoDeleteQueue(bool deleteQueue);

    /**
     * @brief Adds a ProcessQueueTask object to this controller
     * @param t The task to add.
     */
    void addTask(QThread* t);

  protected:

    /**
     * @brief This is the entry point for the task. This is called from the QThread::started signal.
     */
    virtual void run();

  public slots:

    /**
     * @brief Slot that is hooked up to a "finished()" signal from another task. Chaining
     * the signals and slots in this way will ensure that a single new process is begun
     * processing when another task is finished.
     */
    void processTask();

  private:
    QVector<QThread*>  m_Tasks;
    QVector<QThread*>  m_CompletedTasks;
    int m_MaxThreads;
    int m_ThreadCount;
    QMutex              m_Mutex;

    bool            m_AutoDelete;
    ProcessQueueController(const ProcessQueueController&);    // Copy Constructor Not Implemented
    void operator=(const ProcessQueueController&);  // Operator '=' Not Implemented

};

#endif /* PROCESSQUEUECONTROLLER_H_ */
