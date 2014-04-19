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
#ifndef LAYERSDOCKWIDGET_H_
#define LAYERSDOCKWIDGET_H_

#include <QtGui/QDockWidget>

#include "ui_LayersDockWidget.h"
#include "MRCGraphicsView.h"



class LayersDockWidget : public QDockWidget , private Ui::LayersDockWidget
{

    Q_OBJECT;
  public:
    LayersDockWidget(QWidget* parent = 0);
    virtual ~LayersDockWidget();

    void setGraphicsView(MRCGraphicsView* view);

    QCheckBox* getSegmentedImageCheckBox();
    QCheckBox* getOriginalImageCheckBox();
    QSlider*   getOpacitySlider();
    QSpinBox*  getOpacitySpinBox();
    QComboBox* getCompositeTypeComboBox();
    QCheckBox* getUseColorTable();

  public slots:
    void on_segmentedImage_stateChanged(int state);
    void on_originalImage_stateChanged(int state);
    void on_opacitySlider_valueChanged(int value);
    void on_opacitySpinner_valueChanged(int value);
    void on_compositeModeCB_currentIndexChanged();
    void on_useColorTable_stateChanged(int state);

  protected:
    void setupGui();
    void updateDisplayState();


  private:
    MRCGraphicsView* m_GraphicsView;

    LayersDockWidget(const LayersDockWidget&); // Copy Constructor Not Implemented
    void operator=(const LayersDockWidget&); // Operator '=' Not Implemented
};

#endif /* LAYERSDOCKWIDGET_H_ */
