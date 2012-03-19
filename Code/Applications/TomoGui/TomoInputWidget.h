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

 //
 //  This code was written under United States Air Force Contract number
 //                           FA8650-07-D-5800
 //
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef TOMOINPUTWIDGET_H_
#define TOMOINPUTWIDGET_H_


#include <QtCore/QObject>
#include <QtCore/QSettings>
#include <QtGui/QWidget>

#include "ui_TomoInputWidget.h"




/*
 *
 */
class TomoInputWidget : public QWidget, private Ui::TomoInputWidget
{
  Q_OBJECT;

  public:
    TomoInputWidget(QWidget *parent = 0);
    virtual ~TomoInputWidget();

    /**
     * @brief Reads the Preferences from the users pref file
     */
    void readSettings(QSettings &prefs);

    /**
     * @brief Writes the preferences to the users pref file
     */
    void writeSettings(QSettings &prefs);

    void setResolutionMultiple(int x);
    void setIndexLabel(int i);

    double getStopThreshold();
    int getOuterIterations();
    int getInnerIterations();
    double getSigmaX();
    double getMRF();
    int getXYPixelMultiple();
    int getXZPixelMultiple();
    double getUserDefinedOffset();
    bool useDefinedOffset();

  protected slots:


  protected:
    void setupGui();


  private:
    TomoInputWidget(const TomoInputWidget&); // Copy Constructor Not Implemented
    void operator=(const TomoInputWidget&); // Operator '=' Not Implemented
};

#endif /* TOMOINPUTWIDGET_H_ */
