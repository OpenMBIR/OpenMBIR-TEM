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

#include "LayersDockWidget.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LayersDockWidget::LayersDockWidget(QWidget* parent) :
  QDockWidget(parent),
  m_GraphicsView(NULL)
{
  setupUi(this);
  setupGui();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LayersDockWidget::~LayersDockWidget()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::setupGui()
{
  compositeModeCB->blockSignals(true);
  compositeModeCB->insertItem(0, "Exclusion", QVariant(TomoGui_Constants::Exclusion));
  compositeModeCB->insertItem(1, "Difference", QVariant(TomoGui_Constants::Difference));
  compositeModeCB->insertItem(2, "Alpha Blend", QVariant(TomoGui_Constants::Alpha_Blend));
  compositeModeCB->setCurrentIndex(2); // Default to an Alpha Blend;

  compositeModeCB->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::updateDisplayState()
{
  float f = (float)opacitySlider->value() / 100.0;
  m_GraphicsView->setOverlayTransparency(f);

  if (useColorTable->isEnabled() && useColorTable->isChecked())
  {
    //  m_GraphicsView->useCustomColorTable(true);
  }
  else if (useColorTable->isEnabled() && !useColorTable->isChecked())
  {
//   m_GraphicsView->useCustomGrayScaleTable(true);
  }
  else
  {
//    m_GraphicsView->useCustomColorTable(false);
//    m_GraphicsView->useCustomGrayScaleTable(false);
  }

  bool ok = false;
  // Display only the original image only
  if (originalImage->isChecked() && !segmentedImage->isChecked())
  {
    m_GraphicsView->setImageDisplayType(TomoGui_Constants::OriginalImage);
    opacitySlider->setEnabled(false);
    opacitySpinner->setEnabled(false);
    compositeModeCB->setEnabled(false);
    originalImage->setEnabled(false);
  }
  else if ( !originalImage->isChecked() && segmentedImage->isChecked())
  {
    m_GraphicsView->setImageDisplayType(TomoGui_Constants::SegmentedImage);
    opacitySlider->setEnabled(false);
    opacitySpinner->setEnabled(false);
    compositeModeCB->setEnabled(false);
    segmentedImage->setEnabled(false);
  }
  else if (originalImage->isChecked() && segmentedImage->isChecked())
  {
    m_GraphicsView->setImageDisplayType(TomoGui_Constants::CompositedImage);
    TomoGui_Constants::CompositeType cType = static_cast<TomoGui_Constants::CompositeType>(compositeModeCB->itemData(compositeModeCB->currentIndex()).toInt(&ok));
    m_GraphicsView->setCompositeMode(cType);
    compositeModeCB->setEnabled(true);
    if (cType == TomoGui_Constants::Alpha_Blend)
    {
      opacitySlider->setEnabled(true);
      opacitySpinner->setEnabled(true);
    }
    else
    {
      opacitySlider->setEnabled(false);
      opacitySpinner->setEnabled(false);
    }
    originalImage->setEnabled(true);
    segmentedImage->setEnabled(true);
  }
  else
  {
    m_GraphicsView->setImageDisplayType(TomoGui_Constants::UnknownDisplayType);
    opacitySlider->setEnabled(false);
    opacitySpinner->setEnabled(false);
    compositeModeCB->setEnabled(false);
  }
  m_GraphicsView->updateDisplay();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::setGraphicsView(MRCGraphicsView* view)
{
  m_GraphicsView = view;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_segmentedImage_stateChanged(int state)
{
  updateDisplayState();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_originalImage_stateChanged(int state)
{
  updateDisplayState();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_opacitySlider_valueChanged(int value)
{
// std::cout << "on_opacitySlider_valueChanged" << std::endl;
  opacitySpinner->blockSignals(true);
  opacitySpinner->setValue(value);
  opacitySpinner->blockSignals(false);
  float f = (float)value / 100.0;
  m_GraphicsView->setOverlayTransparency(f);
  m_GraphicsView->updateDisplay();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_opacitySpinner_valueChanged(int value)
{
//  std::cout << "on_opacitySpinner_valueChanged" << std::endl;
  opacitySlider->blockSignals(true);
  opacitySlider->setValue(value);
  opacitySlider->blockSignals(false);
  float f = (float)value / 100.0;
  m_GraphicsView->setOverlayTransparency(f);
  m_GraphicsView->updateDisplay();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_compositeModeCB_currentIndexChanged()
{
//  std::cout << "on_compositeModeCB_indexChanged" << std::endl;
  bool ok = false;
  TomoGui_Constants::CompositeType cType = static_cast<TomoGui_Constants::CompositeType>(compositeModeCB->itemData(compositeModeCB->currentIndex()).toInt(&ok));
  m_GraphicsView->setImageDisplayType(TomoGui_Constants::CompositedImage);
  m_GraphicsView->setCompositeMode(cType);
  if (cType == TomoGui_Constants::Alpha_Blend)
  {
    opacitySlider->setEnabled(true);
    opacitySpinner->setEnabled(true);
  }
  else
  {
    opacitySlider->setEnabled(false);
    opacitySpinner->setEnabled(false);
  }
  m_GraphicsView->updateDisplay();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LayersDockWidget::on_useColorTable_stateChanged(int state)
{
  updateDisplayState();
}




QCheckBox* LayersDockWidget::getSegmentedImageCheckBox() { return segmentedImage; }
QCheckBox* LayersDockWidget::getOriginalImageCheckBox() { return originalImage; }
QSlider*   LayersDockWidget::getOpacitySlider() { return opacitySlider; }
QSpinBox*  LayersDockWidget::getOpacitySpinBox() {return opacitySpinner;}
QComboBox* LayersDockWidget::getCompositeTypeComboBox() { return compositeModeCB; }
QCheckBox* LayersDockWidget::getUseColorTable() { return useColorTable; }

