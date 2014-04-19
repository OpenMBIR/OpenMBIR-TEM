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

#include "ProcessQueueDialog.h"

#include <QtCore/QFileInfo>
#include "QProgressLabel.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ProcessQueueDialog::ProcessQueueDialog(QWidget* parent) :
  QWidget(parent)
{
  setupUi(this);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ProcessQueueDialog::~ProcessQueueDialog()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueDialog::clearTable()
{
//  this->processTableWidget->clearContents();
//  this->processTableWidget->setRowCount(0);
  m_TasksMap.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueDialog::addProcess(ProcessQueueTask* task)
{
  verticalLayout->removeItem(verticalSpacer);


  QProgressLabel* progBar = new QProgressLabel(NULL);
  progBar->setRange(0, 100);
  progBar->setAlignment(Qt::AlignBottom);

// QFileInfo fileInfo(task->getInputFilePath());
// progBar->setText(fileInfo.fileName());

  QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
  sizePolicy1.setHorizontalStretch(0);
  sizePolicy1.setVerticalStretch(25);
  sizePolicy1.setHeightForWidth(progBar->sizePolicy().hasHeightForWidth());
  progBar->setSizePolicy(sizePolicy1);

  progBar->setMinimumSize(QSize(1, 25));
  QFont font;
  font.setBold(true);
  font.setWeight(75);
  progBar->setFont(font);

// this->processTableWidget->setRowCount(rowCount + 1);
// this->processTableWidget->setCellWidget(rowCount, 1, progBar);
  verticalLayout->addWidget(progBar);

  verticalLayout->addSpacerItem(verticalSpacer);

  connect(task, SIGNAL(progressValueChanged(int)), progBar, SLOT(setValue(int)));
  connect(task, SIGNAL(taskFinished(QObject*)), this, SLOT(removeRow(QObject*)));
// connect(cancelBtn, SIGNAL(clicked()), task, SLOT(cancel()));
  m_TasksMap[task] = progBar;

// processTableWidget->resizeColumnToContents(0);
// int width = processTableWidget->columnWidth(0);
// width = width + (width * .1);
// processTableWidget->setColumnWidth(0, width);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ProcessQueueDialog::removeRow(QObject* sender)
{
// std::cout << "ProcessQueueDialog::removeRow" << std::endl;
  QWidget* widget = m_TasksMap.value(sender);
  if (NULL != widget)
  {
    verticalLayout->removeWidget(widget);
    m_TasksMap.remove(sender);
    delete widget;
  }

}
