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

#ifndef PROCESSQUEUEDIALOG_H_
#define PROCESSQUEUEDIALOG_H_

#include <ui_ProcessQueueDialog.h>
#include <QtCore/QMap>
#include <QtGui/QWidget>

#include "ProcessQueueTask.h"



/**
 * @class ProcessQueueDialog ProcessQueueDialog.h QtSupport/ProcessQueueDialog.h
 * @brief A QtDialog based class that can display the progress of ProcessQueueTask
 * objects that are running in separate threads as part of a ProcessQueueController.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Jul 26, 2010
 * @version 1.0
 */
class ProcessQueueDialog : public QWidget, private Ui::ProcessQueueDialog
{
    Q_OBJECT
  public:
    ProcessQueueDialog(QWidget *parent = 0);
    virtual ~ProcessQueueDialog();

    /**
     * @brief Adds a ProcessQueueTask to the dialog for monitoring
     * @param task A ProcessQueueTask
     */
    void addProcess(ProcessQueueTask* task);

    /**
     * @brief removes all the ProcessQueueTask Objects from being monitored by
     * the dialog.
     */
    void clearTable();

   signals:
   /**
    * @brief Signal emitted when a ProcessQueueTask is completed.
    */
     void rowComplete(QWidget* widget);

   protected slots:

   /**
    * @brief Will remove a row from the Dialog. Each row represents a single ProcessQueueTask
    * object.
    * @param sender The sender of the message
    */
    void removeRow(QObject* sender);

  private:
    QMap<QObject*, QWidget*> m_TasksMap;

    ProcessQueueDialog(const ProcessQueueDialog&); // Copy Constructor Not Implemented
    void operator=(const ProcessQueueDialog&); // Operator '=' Not Implemented
};

#endif /* PROCESSQUEUEDIALOG_H_ */
