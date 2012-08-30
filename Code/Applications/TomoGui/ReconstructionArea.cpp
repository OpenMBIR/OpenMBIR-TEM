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
#include "ReconstructionArea.h"

#include <iostream>



#include <QtGui/QGraphicsScene>
#include <QtGui/QGraphicsSceneMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QStyleOption>
#include <QtGui/QMenu>
#include <QtGui/QKeyEvent>
#include <QtCore/QRect>


#include "ReconstructionArea.h"
#include "MRCGraphicsView.h"

namespace UIA
{
  const static int Alpha = 155;
}

#define CTRL_POINT_SIZE 10.0f

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionArea::ReconstructionArea(const QPolygonF &polygon, QSize imageSize, QGraphicsItem *parent) :
QGraphicsPolygonItem(polygon, parent),
m_GView(NULL),
m_ImageSize(imageSize),
m_ControlPointMultiplier(1.0f)
{
  setFlag(QGraphicsItem::ItemIsMovable, true);
  setFlag(QGraphicsItem::ItemIsSelectable, true);
  setFlag(QGraphicsItem::ItemIsFocusable, true);
  // setFlag(QGraphicsItem::ItemSendsGeometryChanges, true);
  setAcceptHoverEvents(true);
  m_isResizing = false;
  m_CurrentResizeHandle = ReconstructionArea::NO_CTRL_POINT;
  int maxDim = (imageSize.height() > imageSize.width()) ? imageSize.height() : imageSize.width() ;
  m_CtrlPointSize = maxDim * 0.020; // CTRL_POINT_SIZE * m_ControlPointMultiplier;
  m_GrayLevel = 0 + (255 / 16 * 0);

  m_LineWidth = 1.0;
  m_Visible = true;

  QPoint p = pos().toPoint();
  QRect b = boundingRect().toAlignedRect();
  m_UpperLeft[0] = b.x() + p.x();
  m_UpperLeft[1] = b.y() + p.y();
  m_LowerRight[0] = b.x() + p.x() + b.width();
  m_LowerRight[1] = b.y() + p.y() + b.height();

  QStringList colorNames = QColor::colorNames();
  setBrush(QBrush(QColor(colorNames[21])));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionArea::~ReconstructionArea()
{
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setControlPointMultiplier(float f)
{
    m_ControlPointMultiplier = f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float ReconstructionArea::getControlPointMultiplier()
{
    return m_ControlPointMultiplier;
}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::getUpperLeft( int &x, int &y)
{
  QPointF p = pos();//.toPoint();
  QRect b = boundingRect().toAlignedRect();

  m_UpperLeft[0] = b.x() + p.x();
  m_UpperLeft[1] = b.y() + p.y();

  x = m_UpperLeft[0];
  y = m_UpperLeft[1];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::getLowerRight( int &x, int &y)
{
  QPointF p = pos();//.toPoint();
  QRect b = boundingRect().toAlignedRect();

  m_LowerRight[0] = b.x() + p.x() + b.width();
  m_LowerRight[1] = b.y() + p.y() + b.height();

  x = m_LowerRight[0];
  y = m_LowerRight[1];
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setLineWidth(qreal w)
{
  m_LineWidth = w;
  emit fireReconstructionVOIUpdated(this);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
qreal ReconstructionArea::getLineWidth()
{
  return m_LineWidth;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setVisible(bool visible)
{
  m_Visible = visible;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
  if (m_Visible == false)
  {
    std::cout << "ReconstructionArea::paint Visible=false" << std::endl;
    return;
  }

  painter->setRenderHint(QPainter::Antialiasing, true);

  if (option->state & QStyle::State_HasFocus)
  {
    painter->setPen(QPen(QColor(255, 255, 25, UIA::Alpha), 3.0));
    painter->setBrush(brush());
  }
  else
  {
    painter->setPen(QPen(QColor(255, 25, 25, UIA::Alpha), 3.0));
    painter->setBrush(brush());
  }


  //painter->drawRect(boundingRect());
  float x = boundingRect().x();
  float y = boundingRect().y();
  float w = boundingRect().width();
  float h = boundingRect().height();

  float ctrlPointSize = m_CtrlPointSize; // m_CtrlPointSize * m_ControlPointMultiplier;


//  painter->setPen(QPen(QColor(0, 255, 0, UIA::Alpha), 1.0));
//  QRectF left(x-ctrlPointSize, y, 2*ctrlPointSize, h);
//  QRectF right(x+w-ctrlPointSize, y, 2*ctrlPointSize, h);
//  QRectF top(x, y-ctrlPointSize, w, 2*ctrlPointSize);
//  QRectF bottom(x, y+h-ctrlPointSize, w, 2*ctrlPointSize);
//  painter->drawRect(left);
//  painter->drawRect(right);
//  painter->drawRect(top);
//  painter->drawRect(bottom);


  painter->setPen(QPen(QColor(0, 255, 25, 255), 2.0));
  painter->drawLine(x, y, x+w, y);
  painter->drawLine(x+w, y, x+w, y+h);
  painter->drawLine(x+w, y+h, x, y+h);
  painter->drawLine(x, y+h, x, y);

#if 0
  if (option->state & QStyle::State_Selected)
  {
    painter->setPen(QPen(QColor(0, 255, 25, 255)));
    painter->setBrush( QBrush(QColor(25, 25, 25, UIA::Alpha)));
    //Upper Left
    painter->drawRect((int)x, (int)y, (int)ctrlPointSize, (int)ctrlPointSize);
    //Upper Right
    painter->drawRect((int)x + (int)w - (int)ctrlPointSize, (int)y, (int)ctrlPointSize, (int)ctrlPointSize);
    // Lower Right
    painter->drawRect((int)x + (int)w - (int)ctrlPointSize, (int)y + (int)h - (int)ctrlPointSize, (int)ctrlPointSize, (int)ctrlPointSize);
    // Lower Left
    painter->drawRect((int)x, (int)y + (int)h - (int)ctrlPointSize, (int)ctrlPointSize, (int)ctrlPointSize);
  }
#endif
  painter->setRenderHint(QPainter::Antialiasing, false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::mouseDoubleClickEvent (QGraphicsSceneMouseEvent *event)
{
  if (!isSelected() && scene()) {
      scene()->clearSelection();
      setSelected(true);
  }

//      propertiesSelectedItems(scene());
//    else if (selectedAction == growAction)
//        growSelectedItems(scene());
//    else if (selectedAction == shrinkAction)
//        shrinkSelectedItems(scene());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::contextMenuEvent(QGraphicsSceneContextMenuEvent *event)
{
#if 0
    if (!isSelected() && scene()) {
        scene()->clearSelection();
        setSelected(true);
    }

    QMenu menu;

    QAction *propertiesAction = menu.addAction("Properties");
    menu.addSeparator();
    QAction *delAction = menu.addAction("Delete");


    QAction *selectedAction = menu.exec(event->screenPos());


    if (selectedAction == delAction)
        deleteSelectedItems(scene());
    else if (selectedAction == propertiesAction)
        propertiesSelectedItems(scene());

#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setGraphicsView(MRCGraphicsView* gView)
{
  m_GView = gView;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
 // std::cout << "ReconstructionArea::mousePressEvent" << std::endl;
  static qreal z = 0.0;
  setZValue(z += 1.0);
  m_CurrentResizeHandle = isInResizeArea(event->pos());
  if (event->button() == Qt::LeftButton && m_CurrentResizeHandle != ReconstructionArea::NO_CTRL_POINT)
  {
  //  std::cout << "mousePressEvent m_isResizing = true" << std::endl;
    m_isResizing = true;
  }
  else
  {
    scene()->clearSelection();
    setSelected(true);
    QGraphicsItem::mousePressEvent(event);
    emit fireReconstructionVOISelected(this);
  }

  if (m_GView != NULL)
  {
    m_GView->setAddReconstructionArea(false);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    if (event->button() == Qt::LeftButton && m_isResizing) {
        m_isResizing = false;
    } else {
        QGraphicsItem::mouseReleaseEvent(event);
    }
    emit fireReconstructionVOIUpdated(this);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if (m_isResizing)
  {
    //      std::cout << "mouseMoveEvent m_isResizing = true" << std::endl;
    QPointF lastPos = event->lastScenePos();
    QPointF pos = event->scenePos();
    float deltaX = pos.x() - lastPos.x();
    float deltaY = pos.y() - lastPos.y();
    float x = boundingRect().x();
    float y = boundingRect().y();
    float w = boundingRect().width();
    float h = boundingRect().height();
    //        std::cout << "Delta(): " << deltaX << ", " << deltaY << std::endl;
  //  std::cout << "newRect: " << x << ", " << y << " (" << w << " x " << h << ")" << std::endl;
    QRectF newRect = boundingRect();
    // Move the upper left corner as it is grown
#if 0
    if (m_CurrentResizeHandle == ReconstructionArea::UPPER_LEFT_CTRL_POINT)
    {
      newRect.setX(x + deltaX);
      newRect.setY(y + deltaY);
      newRect.setWidth(w - deltaX);
      newRect.setHeight(h - deltaY);
    }
    else if (m_CurrentResizeHandle == ReconstructionArea::UPPER_RIGHT_CTRL_POINT)
    {
      newRect.setY(y + deltaY);
      newRect.setWidth(w + deltaX);
    }
    else if (m_CurrentResizeHandle == ReconstructionArea::LOWER_LEFT_CTRL_POINT)
    {
      newRect.setX(x + deltaX);
      newRect.setHeight(h + deltaY);
    }
    else if (m_CurrentResizeHandle == ReconstructionArea::LOWER_RIGHT_CTRL_POINT)
    {
      newRect.setWidth(w + deltaX);
      newRect.setHeight(h + deltaY);
    }
    else if (m_CurrentResizeHandle == ReconstructionArea::LEFT_CTRL_POINT)
    {
      newRect.setX(x + deltaX);
      newRect.setWidth(w - deltaX);
    }

    else if (m_CurrentResizeHandle == ReconstructionArea::RIGHT_CTRL_POINT)
    {
      newRect.setWidth(w + deltaX);
    }

    else
#endif
      if (m_CurrentResizeHandle == ReconstructionArea::TOP_CTRL_POINT)
    {
      newRect.setY(y + deltaY);
      newRect.setHeight(h - deltaY);
    }
    else if (m_CurrentResizeHandle == ReconstructionArea::BOTTOM_CTRL_POINT)
    {
      newRect.setHeight(h + deltaY);
    }

  //  std::cout << "A-newRect:" << newRect.x() << ", " << newRect.y() << "   " << newRect.width() << " x " << newRect.height() << std::endl;

    QRectF point = mapRectToItem(m_GView->getImageGraphicsItem(), newRect);
  //  QRect pointInt(point.x(), point.y(), point.width(), point.height());

    // This section sanity checks the rubber band box to make sure it stays in the
    // boundaries of the underlying image
    // Check the x,y corner first and adjust if necessary
    if (point.x() < 0)
    {
      newRect.setX(newRect.x() - point.x());
      point.setX(0.0);
    }
    if (point.y() < 0)
    {
      newRect.setY(newRect.y() - point.y() );
      point.setY(0.0);
    }
  //  std::cout << "point:" << point.x() << ", " << point.y() << "   " << point.width() << " x " << point.height() << std::endl;

    // Now check the width and Height of the image
    //point = mapRectToItem(m_GView->getImageGraphicsItem(), newRect);

    qreal delta = m_ImageSize.width() - (point.x() + point.width());
    if (delta < 0.0 )
    {
  //    std::cout << "delta X:" << delta << std::endl;
      newRect.setWidth(newRect.width() + delta);
    }
    delta = m_ImageSize.height() - (point.y() + point.height());
    if (delta < 0.0 )
    {
  //    std::cout << "delta Y:" << delta << std::endl;
      newRect.setHeight(newRect.height() + delta);
    }

  //  std::cout << "B-newRect:" << newRect.x() << ", " << newRect.y() << "   " << newRect.width() << " x " << newRect.height() << std::endl;
    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
  }
  else
  {
   // QGraphicsItem::mouseMoveEvent(event);
  }
  emit fireReconstructionVOIUpdated(this);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::updateGeometry(int xmin, int ymin, int xmax, int ymax)
{
    std::cout << "xmin\tymin\txmax\tymax" << std::endl;
    //std::cout << xmin << "\t" << ymin << "\t" << xmax << "\t" << ymax << std::endl;
    // Need to flip the Y's because we display our image origin in the lower left
    qint32 yTemp = ymin;
    ymin = m_ImageSize.height() - ymax;
    ymax = m_ImageSize.height() - yTemp;
    std::cout << xmin << "\t" << ymin << "\t" << xmax << "\t" << ymax << std::endl;

    QRectF newRect(xmin, ymin, (xmax - xmin), (ymax - ymin));

    getUpperLeft(xmin, ymin);
    getLowerRight(xmax, ymax);
    std::cout << xmin << "\t" << ymin << "\t" << xmax << "\t" << ymax << std::endl;

    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::updateWidth(float percentWidth)
{
  float remWidth = m_ImageSize.width() * percentWidth/2.0;
  float midWidth = m_ImageSize.width()/2.0f;

  qint32 xmin = 0;
  qint32 xmax = 0;
  qint32 ymin = 0;
  qint32 ymax = 0;
  getUpperLeft(xmin, ymin);
  getLowerRight(xmax, ymax);

  float x_min = midWidth - remWidth;
  float x_max = midWidth + remWidth;

  QRectF newRect(x_min, ymin, (x_max - x_min), (ymax - ymin));
  prepareGeometryChange();
  setPolygon(QPolygonF(newRect));

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setXMin(const QString &str)
{
  qint32 xmin = 0;
  qint32 xmax = 0;
  qint32 ymin = 0;
  qint32 ymax = 0;
  getUpperLeft(xmin, ymin);
  getLowerRight(xmax, ymax);

  bool ok = false;
  qint32 x_min = str.toInt(&ok);
 // if (x_min < xmax)
  {
    QRectF newRect(x_min, ymin, (xmax - x_min), (ymax - ymin));
    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setYMin(const QString &str)
{
  qint32 xmin = 0;
  qint32 xmax = 0;
  qint32 ymin = 0;
  qint32 ymax = 0;
  getUpperLeft(xmin, ymin);
  getLowerRight(xmax, ymax);

  bool ok = false;
  qint32 y_min = m_ImageSize.height() - str.toInt(&ok);

 // if (y_min < ymax)
  {
    QRectF newRect(xmin, y_min, (xmax - xmin), (ymax - y_min));
    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setXMax(const QString &str)
{
  qint32 xmin = 0;
  qint32 xmax = 0;
  qint32 ymin = 0;
  qint32 ymax = 0;
  getUpperLeft(xmin, ymin);
  getLowerRight(xmax, ymax);
  bool ok = false;
  qint32 x_max = str.toInt(&ok);
 // if (x_max > xmin)
  {
    QRectF newRect(xmin, ymin, (x_max - xmin), (ymax - ymin));
    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setYMax(const QString &str)
{
  qint32 xmin = 0;
  qint32 xmax = 0;
  qint32 ymin = 0;
  qint32 ymax = 0;
  getUpperLeft(xmin, ymin);
  getLowerRight(xmax, ymax);

  bool ok = false;
  qint32 y_max = m_ImageSize.height() - str.toInt(&ok);
 // if (y_max > ymin)
  {
    QRectF newRect(xmin, ymin, (xmax - xmin), (y_max - ymin));
    prepareGeometryChange();
    setPolygon(QPolygonF(newRect));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::hoverMoveEvent(QGraphicsSceneHoverEvent *event)
{
  ReconstructionArea::CTRL_POINTS pt = isInResizeArea(event->pos());
  if ( (m_isResizing || pt != ReconstructionArea::NO_CTRL_POINT) && isSelected()  )
  {

    if (pt == UPPER_LEFT_CTRL_POINT || pt == LOWER_RIGHT_CTRL_POINT)
    {
      setCursor(Qt::SizeFDiagCursor);
    }
    else if (pt == UPPER_RIGHT_CTRL_POINT || pt == LOWER_LEFT_CTRL_POINT)
    {
      setCursor(Qt::SizeBDiagCursor);
    }
    else if (pt == TOP_CTRL_POINT || pt == BOTTOM_CTRL_POINT)
    {
      setCursor(Qt::SizeVerCursor);
    }
    else if (pt == LEFT_CTRL_POINT || pt == RIGHT_CTRL_POINT)
    {
      setCursor(Qt::SizeHorCursor);
    }
  }
  else
  {
    setCursor(Qt::ArrowCursor);
  }
  QGraphicsItem::hoverMoveEvent(event);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Delete:
        break;
    case Qt::Key_Insert:
        break;
    case Qt::Key_Plus:
        break;
    case Qt::Key_Minus:
        break;
    default:
        QGraphicsItem::keyPressEvent(event);
        break;
    }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionArea::CTRL_POINTS ReconstructionArea::isInResizeArea(const QPointF &pos)
{
  float x = boundingRect().x();
  float y = boundingRect().y();
  float w = boundingRect().width();
  float h = boundingRect().height();

  float ctrlPointSize = m_CtrlPointSize; // m_CtrlPointSize * m_ControlPointMultiplier;
  QRectF upLeft(x, y, ctrlPointSize, ctrlPointSize);
  QRectF upRight(x + w - ctrlPointSize, y, ctrlPointSize, ctrlPointSize);
  QRectF lowRight(x + w - ctrlPointSize, y + h - ctrlPointSize, ctrlPointSize, ctrlPointSize);
  QRectF lowLeft(x, y + h - ctrlPointSize, ctrlPointSize, ctrlPointSize);

  QRectF left(x-ctrlPointSize, y, 2*ctrlPointSize, h);
  QRectF right(x+w-5, y, 2*ctrlPointSize, h);
  QRectF top(x, y-ctrlPointSize, w, 2*ctrlPointSize);
  QRectF bottom(x, y+h-ctrlPointSize, w, 2*ctrlPointSize);
#if 0
  if (upLeft.contains(pos))
  {
   // std::cout << "UPPER_LEFT_CTRL_POINT" << std::endl;
    return UPPER_LEFT_CTRL_POINT;
  }
  if (upRight.contains(pos))   {
   // std::cout << "UPPER_RIGHT_CTRL_POINT" << std::endl;
    return UPPER_RIGHT_CTRL_POINT;
  }
  if (lowRight.contains(pos))   {
  //  std::cout << "LOWER_RIGHT_CTRL_POINT" << std::endl;
    return LOWER_RIGHT_CTRL_POINT;
  }
  if (lowLeft.contains(pos))   {
   // std::cout << "LOWER_LEFT_CTRL_POINT" << std::endl;
    return LOWER_LEFT_CTRL_POINT;
  }
  if (left.contains(pos)){
  //  std::cout << "LEFT_CTRL_POINT" << std::endl;
    return LEFT_CTRL_POINT;
  }
  if (right.contains(pos)) {
  //  std::cout << "RIGHT_CTRL_POINT" << std::endl;
    return RIGHT_CTRL_POINT;
  }
#endif
  if (top.contains(pos)) {
  //  std::cout << "TOP_CTRL_POINT" << std::endl;
    return TOP_CTRL_POINT;
  }

  if (bottom.contains(pos)) {
  //  std::cout << "BOTTOM_CTRL_POINT" << std::endl;
    return BOTTOM_CTRL_POINT;
  }

  return NO_CTRL_POINT;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ReconstructionArea::type() const
{
    return Type;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReconstructionArea::setColor(QColor color)
{
  setBrush(QBrush(color));
  emit fireReconstructionVOIUpdated(this);
}

