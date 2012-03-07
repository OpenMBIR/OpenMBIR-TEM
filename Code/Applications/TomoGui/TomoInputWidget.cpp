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
#include "TomoInputWidget.h"

#include <iostream>

#define READ_STRING_SETTING(prefs, var, emptyValue)\
  var->setText( prefs.value(#var).toString() );\
  if (var->text().isEmpty() == true) { var->setText(emptyValue); }


#define READ_SETTING(prefs, var, ok, temp, default, type)\
  ok = false;\
  temp = prefs.value(#var).to##type(&ok);\
  if (false == ok) {temp = default;}\
  var->setValue(temp);

#define READ_VALUE(prefs, var, ok, temp, default, type)\
  ok = false;\
  temp = prefs.value(#var).to##type(&ok);\
  if (false == ok) {temp = default;}\
  var = temp;

#define WRITE_STRING_SETTING(prefs, var)\
  prefs.setValue(#var , this->var->text());

#define WRITE_SETTING(prefs, var)\
  prefs.setValue(#var, this->var->value());

#define READ_BOOL_SETTING(prefs, var, emptyValue)\
  { QString s = prefs.value(#var).toString();\
  if (s.isEmpty() == false) {\
    bool bb = prefs.value(#var).toBool();\
  var->setChecked(bb); } else { var->setChecked(emptyValue); } }

#define WRITE_BOOL_SETTING(prefs, var, b)\
    prefs.setValue(#var, (b) );

#define WRITE_CHECKBOX_SETTING(prefs, var)\
    prefs.setValue(#var, var->isChecked() );

#define WRITE_VALUE(prefs, var)\
    prefs.setValue(#var, var);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoInputWidget::TomoInputWidget(QWidget *parent) :
QWidget(parent)
{
  setupUi(this);
  setupGui();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoInputWidget::~TomoInputWidget()
{

}

// -----------------------------------------------------------------------------
//  Read the prefs from the local storage file
// -----------------------------------------------------------------------------
void TomoInputWidget::readSettings(QSettings &prefs)
{
  QString val;
  bool ok;
  qint32 i;
  double d;
  prefs.beginGroup(interationIndex->text());
  std::cout << "Reading Settings for Index " << interationIndex->text().toStdString() << std::endl;


  READ_STRING_SETTING(prefs, interationIndex, "");
  READ_STRING_SETTING(prefs, stopThreshold, "0.01")
  READ_STRING_SETTING(prefs, innerIterations, "10")
  READ_STRING_SETTING(prefs, outerIterations, "1")
  READ_STRING_SETTING(prefs, sigmaX, "0.00004")
  READ_SETTING(prefs, mrf, ok, d, 0.00001, Double)
  READ_SETTING(prefs, xyPixelMultiple, ok, i, 1, Int)
  READ_SETTING(prefs, xzPixelMultiple, ok, i, 1, Int)
  READ_BOOL_SETTING(prefs, useDefaultOffset, false)
  READ_STRING_SETTING(prefs, defaultOffset, "0")
  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//  Write our prefs to file
// -----------------------------------------------------------------------------
void TomoInputWidget::writeSettings(QSettings &prefs)
{
  QString val;
  bool ok;
//  qint32 i;
//  double d;

  prefs.beginGroup(interationIndex->text());
  //std::cout << "Writing Settings for Index " << interationIndex->text().toStdString() << std::endl;


  WRITE_STRING_SETTING(prefs, interationIndex);
  WRITE_STRING_SETTING(prefs, stopThreshold)
  WRITE_STRING_SETTING(prefs, innerIterations)
  WRITE_STRING_SETTING(prefs, outerIterations)
  WRITE_STRING_SETTING(prefs, sigmaX)
  WRITE_SETTING(prefs, mrf)
  WRITE_SETTING(prefs, xyPixelMultiple)
  WRITE_SETTING(prefs, xzPixelMultiple)
  WRITE_BOOL_SETTING(prefs, useDefaultOffset, useDefaultOffset->isChecked())
  WRITE_STRING_SETTING(prefs, defaultOffset)
  prefs.endGroup();
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoInputWidget::setupGui()
{
  QDoubleValidator* dVal = new QDoubleValidator(this);
  stopThreshold->setValidator(dVal);

  QIntValidator* outVal = new QIntValidator(this);
  outerIterations->setValidator(outVal);

  QIntValidator* inVal = new QIntValidator(this);
  innerIterations->setValidator(inVal);

  QDoubleValidator* sigVal = new QDoubleValidator(this);
  sigmaX->setValidator(sigVal);

  QDoubleValidator* offsetVal = new QDoubleValidator(this);
  defaultOffset->setValidator(offsetVal);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoInputWidget::setIndexLabel(int i)
{
  interationIndex->setText(QString::number(i, 10));
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoInputWidget::setResolutionMultiple(int x)
{
  xyPixelMultiple->setValue(x);
  xzPixelMultiple->setValue(x);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double TomoInputWidget::getStopThreshold()
{
  bool ok = false;
  return stopThreshold->text().toDouble(&ok);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TomoInputWidget::getOuterIterations()
{
  bool ok = false;
  return outerIterations->text().toInt(&ok);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TomoInputWidget::getInnerIterations()
{
  bool ok = false;
  return innerIterations->text().toInt(&ok);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double TomoInputWidget::getSigmaX()
{
  bool ok = false;
  return sigmaX->text().toDouble(&ok);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double TomoInputWidget::getMRF()
{
  return mrf->value();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TomoInputWidget::getXYPixelMultiple()
{
  return xyPixelMultiple->value();
}

// ----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TomoInputWidget::getXZPixelMultiple()
{
  return xzPixelMultiple->value();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double TomoInputWidget::getUserDefinedOffset()
{
  bool ok = false;
  return defaultOffset->text().toDouble(&ok);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TomoInputWidget::useDefinedOffset()
{
  return useDefaultOffset->isChecked();
}
