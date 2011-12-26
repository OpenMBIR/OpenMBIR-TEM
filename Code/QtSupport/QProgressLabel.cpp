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
#include "QProgressLabel.h"
#include <QtGui/QPen>
#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QBrush>
#include <QtGui/QLinearGradient>

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QProgressLabel::QProgressLabel(QWidget* parent) :
  QLabel(parent),
  m_Color(Qt::blue),
  m_Highlight(Qt::lightGray),
  m_Shadow(Qt::blue),
  m_Opacity(0.70),
  m_Roundness(0),
  m_Maximum(100),
  m_Minimum(0),
  m_Value(0)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QProgressLabel::~QProgressLabel()
{
 // std::cout << "~QProgressLabel" << std::endl;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QProgressLabel::paintEvent(QPaintEvent* event)
{
  Q_UNUSED(event);
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);

  //test for state changes
  QColor button_color = m_Color;

  QRect button_rect = this->geometry();
  QSize size = this->size();
  int width = this->size().width() * static_cast<qreal>(m_Value - m_Minimum) / static_cast<qreal>(m_Maximum - m_Minimum);



  if (m_Value > 0) {
    //outline
    painter.setPen(QPen(QBrush(Qt::darkBlue), 2.0));
    QPainterPath outline;
    outline.addRoundedRect(0, 0, width, button_rect.height(), m_Roundness, m_Roundness);
    painter.setOpacity(m_Opacity);
    painter.drawPath(outline);

    //gradient
    QLinearGradient gradient(0, 0, 0, button_rect.height());
    gradient.setSpread(QGradient::ReflectSpread);
    gradient.setColorAt(0.0, button_color);
    gradient.setColorAt(0.4, m_Shadow);
    gradient.setColorAt(0.6, m_Shadow);
    gradient.setColorAt(1.0, button_color);

    QBrush brush(gradient);
    painter.setBrush(brush);
    painter.setPen(QPen(QBrush(button_color), 2.0));

    //main button
    QPainterPath painter_path;
    painter_path.addRoundedRect(1, 1, width - 1, button_rect.height() - 1, m_Roundness, m_Roundness);
    painter.setClipPath(painter_path);

    painter.setOpacity(m_Opacity);
    painter.drawRoundedRect(1, 1, width - 1, button_rect.height() - 1, m_Roundness, m_Roundness);

    //glass highlight
    painter.setBrush(QBrush(Qt::white));
    painter.setPen(QPen(QBrush(Qt::white), 0.01));
    painter.setOpacity(0.30);
    painter.drawRect(1, 1, width - 1, (button_rect.height() / 2) - 1);
  }
  //text
  QString text = this->text();
  if (!text.isNull())
  {
    QPainterPath painter_path;
    painter_path.addRoundedRect(1, 1, button_rect.width(), button_rect.height() - 1, m_Roundness, m_Roundness);
    painter.setClipPath(painter_path);

    QFont font = this->font();
    painter.setFont(font);
    painter.setPen(Qt::black);
    painter.setOpacity(1.0);
    painter.drawText(0, 0, button_rect.width(), button_rect.height(), Qt::AlignCenter, text);
  }

}
