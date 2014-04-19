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

#ifndef QPROGRESSLABEL_H_
#define QPROGRESSLABEL_H_

#include <QtGui/QLabel>
#include <QtGui/QColor>

/**
* @brief Creates a "setter" method to set the property.
*/
#define QT_SET_PROPERTY(type, prpty) \
  void set##prpty(type value) { this->m_##prpty = value; }

#define QT_GET_PROPERTY(type, prpty, var) \
  type var() { return m_##prpty; }

#define QT_INSTANCE_PROPERTY(type, prpty, var)\
  private:\
  type   m_##prpty;\
  public:\
  QT_SET_PROPERTY(type, prpty)\
  QT_GET_PROPERTY(type, prpty, var)

#define QT_SLOT_PROPERTY(type, prpty, var)\
  private:\
  type   m_##prpty;\
  public:\
  QT_GET_PROPERTY(type, prpty, var)\
   
#include <iostream>
/*
 *
 */
class QProgressLabel : public QLabel
{
    Q_OBJECT;

  public:
    QProgressLabel(QWidget* parent = 0 );
    virtual ~QProgressLabel();



    QT_INSTANCE_PROPERTY(QColor, Color, color);
    QT_INSTANCE_PROPERTY(QColor, Highlight, highlight);
    QT_INSTANCE_PROPERTY(QColor, Shadow, shadow);
    QT_INSTANCE_PROPERTY(qreal, Opacity, opacity);
    QT_INSTANCE_PROPERTY(int, Roundness, roundness);

    QT_SLOT_PROPERTY(int, Maximum, maximum);
    QT_SLOT_PROPERTY(int, Minimum, minimum);
    QT_SLOT_PROPERTY(int, Value, value);

    void setRange(int min, int max) { m_Minimum = min; m_Maximum = max; }


  public slots:
    void setMaximum(int m) { m_Maximum = m; }
    void setMinimum(int m) { m_Minimum = m; }
    void setValue(int v)
    {
      //  std::cout << "Setting Value: " << v << std::endl;
      m_Value = v;
      repaint();
      emit valueChanged(m_Value);
    }

  signals:

    void  valueChanged ( int value );

  protected:

    void paintEvent(QPaintEvent* event);


  private:
    QProgressLabel(const QProgressLabel&); // Copy Constructor Not Implemented
    void operator=(const QProgressLabel&); // Operator '=' Not Implemented

};

#endif /* QPROGRESSLABEL_H_ */
