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

#include "ProcessQueueController.h"
#include <QtCore/QMutexLocker>
#include <iostream>

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ProcessQueueController::ProcessQueueController(QObject* parent) :
  QThread(parent),
  m_AutoDelete(true)
{
  m_MaxThreads = QThread::idealThreadCount();
  m_ThreadCount = 1; // We need this to be 1 because the first time we will decrement the value
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ProcessQueueController::~ProcessQueueController()
{
  for (int i = 0; i < m_CompletedTasks.count(); ++i)
  {
    QThread* t = m_CompletedTasks.at(i);
    t->deleteLater(); // Clean up the memory
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueController::run()
{
  exec(); // <== Starts the Event Loop for this thread. Will BLOCK here until the Quit Slot is called
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueController::setAutoDeleteQueue(bool deleteQueue)
{
  this->m_AutoDelete = deleteQueue;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueController::addTask(QThread* t)
{
  if (NULL == t) { return; }
  if (m_Tasks.contains( t) == false)
  {
    m_Tasks.push_back(t);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueController::processTask()
{
  QMutexLocker lock(&m_Mutex);

  --m_ThreadCount; //decrement the value as this is called when another thread gets finished
//  std::cout << "ProcessQueueController::processTask() Start m_ThreadCount: " << m_ThreadCount << std::endl;
  while (m_ThreadCount < m_MaxThreads)
  {
    if (m_Tasks.count() > 0)
    {
      QThread* task = m_Tasks.front();
      m_Tasks.pop_front(); // Remove the thread from the QVector
      m_CompletedTasks.push_back(task);
      connect(task, SIGNAL(finished()), this, SLOT(processTask()));
      task->start();
      m_ThreadCount++;
    }
    else
    {
      break;
    }
  }
  if (m_ThreadCount == 0)
  {
    //  std::cout << "Last Thread has finished." << std::endl;
    quit(); // Tells the event loop to exit, which will then exit the thread.
  }
}
