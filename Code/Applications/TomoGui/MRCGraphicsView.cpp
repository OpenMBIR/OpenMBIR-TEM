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


#include "MRCGraphicsView.h"

#include <QtCore/QFileInfo>
#include <QtCore/QUrl>
#include <QtGui/QWidget>
#include <QDragEnterEvent>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QtGui/QPixmap>
#include <QtGui/QGraphicsPolygonItem>

#include "TomoGui.h"
#include "ReconstructionArea.h"

namespace UIA
{
  const static int Alpha = 155;
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCGraphicsView::MRCGraphicsView(QWidget *parent)
: QGraphicsView(parent),
  m_ImageGraphicsItem(NULL),
  m_DisableVOISelection(false),
  m_ReconstructionArea(NULL)
{
  setAcceptDrops(true);
  setDragMode(RubberBandDrag);
  m_AddReconstructionArea = true;

  m_ZoomFactors[0] = 0.1f;
  m_ZoomFactors[1] = 0.25f;
  m_ZoomFactors[2] = 0.5f;
  m_ZoomFactors[3] = 1.0f;
  m_ZoomFactors[4] = 1.250f;
  m_ZoomFactors[5] = 1.500f;
  m_ZoomFactors[6] = 2.000f;
  m_ZoomFactors[7] = 4.000f;
  m_ZoomFactors[8] = 6.000f;
  m_ZoomFactors[9] = -1.0f;
  m_MainGui = NULL;
  m_RubberBand = NULL;
  m_ImageDisplayType = TomoGui_Constants::OriginalImage;
  m_composition_mode = QPainter::CompositionMode_SourceOver;
  m_OverlayTransparency = 1.0f; // Fully opaque
//  setBackgroundBrush(QPixmap(":/background4.png"));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setOverlayTransparency(float f)
{
  m_OverlayTransparency = f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::addReconstructionArea(bool b)
{
  m_AddReconstructionArea = b;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::disableVOISelection(bool b)
{
  m_DisableVOISelection = b;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::zoomIn()
{
  scale(1.1, 1.1);
//  if (NULL != m_ReconstructionArea)
//  {
//    m_ReconstructionArea->setControlPointMultiplier(m_ReconstructionArea->getControlPointMultiplier() * 1.1f);
//  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::zoomOut()
{
  scale(1.0 / 1.1, 1.0 / 1.1);
//  if (NULL != m_ReconstructionArea)
//  {
//    m_ReconstructionArea->setControlPointMultiplier(m_ReconstructionArea->getControlPointMultiplier() / 1.1f);
//  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setZoomIndex(int index)
{
  if (index == 9)
  {
    QRectF r = scene()->sceneRect();
    fitInView(r, Qt::KeepAspectRatio);
  }
  else
  {
    QTransform transform;
    transform.scale(m_ZoomFactors[index], m_ZoomFactors[index]);
    setTransform(transform);
//    if (NULL != m_ReconstructionArea)
//    {
//        m_ReconstructionArea->setControlPointMultiplier(m_ReconstructionArea->getControlPointMultiplier() * m_ZoomFactors[index]);
//    }
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::dragEnterEvent(QDragEnterEvent *event)
{
 // qWarning("QFSDroppableGraphicsView::dragEnterEvent(QDragEnterEvent *event)");
  // accept just text/uri-list mime format
  if (event->mimeData()->hasFormat("text/uri-list"))
  {
    event->acceptProposedAction();
  }
  this->setStyleSheet("border: 1px solid green;");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::dragLeaveEvent(QDragLeaveEvent *event)
{
//  qWarning("QFSDroppableGraphicsView::dragLeaveEvent(QDragLeaveEvent *event)");
  this->setStyleSheet("");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::dropEvent(QDropEvent *event)
{
  this->setStyleSheet("");
//  qWarning("QFSDroppableGraphicsView::dropEvent(QDropEvent *event)");
  QList<QUrl> urlList;
  QString fName;
  QFileInfo info;
#if 0
  if (event->mimeData()->hasUrls())
  {
    urlList = event->mimeData()->urls(); // returns list of QUrls
    // if just text was dropped, urlList is empty (size == 0)

    if ( urlList.size() > 0) // if at least one QUrl is present in list
    {
      fName = urlList[0].toLocalFile(); // convert first QUrl to local path
      info.setFile( fName ); // information about file
      QString ext = info.suffix().toLower();
      if (ext.compare("tif") == 0
          || ext.compare("tiff") == 0
          || ext.compare("jpg") == 0
          || ext.compare("jpeg") == 0
          || ext.compare("png") == 0
          || ext.compare("bmp") == 0)
      {
        m_MainGui->on_inputMRCFilePath_textChanged(fName);
      }
    }
  }
#endif
  event->acceptProposedAction();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::updateDisplay()
{
  QPainter painter;
  QSize pSize(0, 0);
  if (m_BaseImage.isNull() == false)
  {
   pSize = m_BaseImage.size();
  }
  else
  {
    return;
  }

  QImage paintImage(pSize, QImage::Format_ARGB32_Premultiplied);
  QPoint point(0, 0);
  painter.begin(&paintImage);
  painter.setPen(Qt::NoPen);
  if (m_ImageDisplayType == TomoGui_Constants::OriginalImage)
  {
    painter.drawImage(point, m_BaseImage);
  }

  painter.end();

  if (paintImage.isNull() == true)
  {
    return;
  }
  QGraphicsPixmapItem *pixItem = qgraphicsitem_cast<QGraphicsPixmapItem*> (m_ImageGraphicsItem);
  pixItem->setPixmap(QPixmap::fromImage(paintImage));

  this->update();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setImageDisplayType(TomoGui_Constants::ImageDisplayType displayType)
{
  m_ImageDisplayType = displayType;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::clearContent()
{
  QGraphicsScene* gScene = scene();

  if (gScene != NULL)
  {

    gScene->removeItem(m_ImageGraphicsItem);
    if(NULL != m_ReconstructionArea)
    {
      m_ReconstructionArea->setParentItem(NULL);
    }
    delete m_ImageGraphicsItem; // Will delete its children
    m_ImageGraphicsItem = NULL;
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::loadBaseImageFile(QImage image)
{
  m_BaseImage = image;
  if(m_BaseImage.isNull() == true)
  {
 //   std::cout << "MRCGraphicsView::loadBaseImageFile() - Input Image was NULL" << std::endl;
    return;
  }
  QSize pSize(0, 0);

  m_BaseImage = m_BaseImage.convertToFormat(QImage::Format_ARGB32_Premultiplied);

  QGraphicsScene* gScene = scene();
  if(gScene == NULL)
  {
    gScene = new QGraphicsScene(this);
    setScene(gScene);
  }
  else
  {
    gScene->removeItem(m_ImageGraphicsItem);
    if(NULL != m_ReconstructionArea)
    {
      m_ReconstructionArea->setParentItem(NULL);
    }
    delete m_ImageGraphicsItem; // Will delete its children
    m_ImageGraphicsItem = NULL;
  }

  m_ImageGraphicsItem = gScene->addPixmap(QPixmap::fromImage(m_BaseImage));
  if(NULL != m_ReconstructionArea)
  {
    m_ReconstructionArea->setParentItem(m_ImageGraphicsItem);
  }

  m_ImageGraphicsItem->setAcceptDrops(true);
  m_ImageGraphicsItem->setZValue(-1);
  QRectF rect = m_ImageGraphicsItem->boundingRect();
  gScene->setSceneRect(rect);

  // Line Color
  m_XZLine.setPen(QPen(QColor(255, 255, 0, UIA::Alpha)));
  // Fill Color
  m_XZLine.setBrush(QBrush(QColor(255, 255, 0, UIA::Alpha)));
  m_XZLine.setZValue(1);
  m_XZLine.setCacheMode(QGraphicsItem::DeviceCoordinateCache);
  m_XZLine.setVisible(m_XZLine.isVisible());
  if(scene()->items().contains(&m_XZLine) == false)
  {
    scene()->addItem(&m_XZLine);
  }

  this->updateDisplay();
  emit fireBaseMRCFileLoaded();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QGraphicsItem* MRCGraphicsView::getImageGraphicsItem()
{
  return m_ImageGraphicsItem;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage MRCGraphicsView::getBaseImage()
{
  return m_BaseImage;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setAddReconstructionArea(bool b)
{
  m_AddReconstructionArea = b;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::mousePressEvent(QMouseEvent *event)
{
  // std::cout << "TomoGuiGraphicsView::mousePressEvent accepted:" << (int)(event->isAccepted()) << std::endl;
  m_AddReconstructionArea = true;
  QGraphicsView::mousePressEvent(event);
  if (m_DisableVOISelection == true)
  {
    return;
  }

  // This next code makes sure the rubber band is within the boundaries of the image
  // and if it is not then it forcefully sets it there. It only checks if we
  // are to the left or above the image. If we are right or below the image
  // the mouseRelease event will take care of that sanity check.
  m_MouseClickOrigin = event->pos();
  QRect box = QRect(m_MouseClickOrigin, event->pos()).normalized();
  QPolygonF sceneBox = mapToScene(box);
  QRectF rectf = sceneBox.boundingRect();

  if (rectf.x() < 0 )
  {
    m_MouseClickOrigin.setX( m_MouseClickOrigin.x() + -rectf.x() );
  }
  if (rectf.y() < 0 )
  {
    m_MouseClickOrigin.setY( m_MouseClickOrigin.y() + -rectf.y() );
  }

  // std::cout << "    event->accepted() == false" << std::endl;
  if(m_AddReconstructionArea == true)
  {
 //   m_MouseClickOrigin = event->pos();
    if(!m_RubberBand) m_RubberBand = new QRubberBand(QRubberBand::Rectangle, this);
    m_RubberBand->setGeometry(QRect(m_MouseClickOrigin, QSize()));
    m_RubberBand->show();
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::mouseMoveEvent(QMouseEvent *event)
{
  if(m_AddReconstructionArea == true && m_RubberBand != NULL)
  {
    m_RubberBand->setGeometry(QRect(m_MouseClickOrigin, event->pos()).normalized());
  }
  else
  {
    QGraphicsView::mouseMoveEvent(event);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::updateXZLine(float percentWidth)
{
    float remWidth = m_BaseImage.size().width() * percentWidth/2.0;

    float midWidth = m_BaseImage.size().width()/2.0f;

    QLineF currentLine = getXZPlane();
    float xStart = midWidth - remWidth;
    float xEnd = midWidth + remWidth;

    QPointF p0(xStart, currentLine.y1());
    QPointF p1(xEnd, currentLine.y1());


    QVector<QPointF> line;
    line.push_back(p0);
    line.push_back(p1);
    QPolygonF polygon(line);
    m_XZLine.setPolygon(polygon);
    m_XZLine.setVisible(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::mouseReleaseEvent(QMouseEvent *event)
{

  if (event->modifiers() == Qt::ShiftModifier)
  {
    m_RubberBand->hide();
    QRectF sr = sceneRect();

    QPointF mappedPoint = mapToScene(QPoint(0, m_MouseClickOrigin.y()));

    if (mappedPoint.y() > m_BaseImage.size().height())
    {
      mappedPoint.setY(m_BaseImage.size().height());
    }

    float percentWidth = 0.9; // Start with a 90% width for single slice reconstructions
    float remWidth = m_BaseImage.size().width() * percentWidth/2.0;
    float midWidth = m_BaseImage.size().width()/2.0f;
    float xStart = midWidth - remWidth;
    float xEnd = midWidth + remWidth;

    QPointF p0(xStart, mappedPoint.y());
    QPointF p1(xEnd, mappedPoint.y());

  //  std::cout << "SceneRect: " << sr.x() << ", " << sr.y() << "  " << sr.width() << ", " << sr.height() << std::endl;
    QVector<QPointF> line;
    line.push_back(p0);
    line.push_back(p1);
    QPolygonF polygon(line);

    m_XZLine.setPolygon(polygon);
    m_XZLine.setVisible(true);

    int y = m_BaseImage.size().height() - getXZPlane().y1();


    emit fireSingleSliceSelected(y);
  }
  else if (m_AddReconstructionArea == true && NULL != m_RubberBand)
  {
    m_RubberBand->hide();
    QPoint endPoint = event->pos();

    QRect box = QRect(m_MouseClickOrigin, event->pos()).normalized();

    QPolygonF sceneBox = mapToScene(box);
    QRectF boxf = sceneBox.boundingRect();
    //make sure we are within the upper left corner of the image
    if (boxf.left() < 0) { boxf.setLeft(0); }
    if (boxf.top() < 0) { boxf.setTop(0); }


    if (boxf.right() > m_BaseImage.size().width())
    {
      boxf.setRight(m_BaseImage.size().width());
    }
    if (boxf.bottom() > m_BaseImage.size().height())
    {
      boxf.setBottom(m_BaseImage.size().height());
    }


    createNewReconstructionArea(boxf);
    m_AddReconstructionArea = false;
  }
  else
  {
    QGraphicsView::mouseReleaseEvent(event);
    if (scene())
    {
      QList<QGraphicsItem *> selected;
      selected = scene()->selectedItems();
      if (selected.count() == 0)
      {
        emit fireReconstructionVOILostFocus();
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QLineF MRCGraphicsView::getXZPlane()
{
 QPointF p0 = m_XZLine.polygon().at(0);
 QPointF p1 = m_XZLine.polygon().at(1);

 QLineF line(p0, p1);
 return line;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::createNewReconstructionArea(const QRectF brect)
{

  ReconstructionArea* userInitArea = new ReconstructionArea(brect, m_BaseImage.size());
  userInitArea->setGraphicsView(this);
  // Line Color
  userInitArea->setPen(QPen(QColor(225, 225, 225, UIA::Alpha)));
  // Fill Color
  userInitArea->setBrush(QBrush(QColor(28, 28, 200, UIA::Alpha)));
  userInitArea->setZValue(1);
  userInitArea->setCacheMode(QGraphicsItem::DeviceCoordinateCache);

  addNewReconstructionArea(userInitArea);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::addNewReconstructionArea(ReconstructionArea* reconVOI)
{
  // Set the Parent Item
  reconVOI->setParentItem(m_ImageGraphicsItem);
  // Add it to the vector of UserInitAreas
  if (m_ReconstructionArea != NULL)
  {
    m_ReconstructionArea->deleteLater();
  }
  m_ReconstructionArea = reconVOI;

  emit fireReconstructionVOIAdded(reconVOI);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReconstructionArea* MRCGraphicsView::reconstructionArea()
{
  return m_ReconstructionArea;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setWidget(QWidget* gui)
{
  m_MainGui = gui;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCGraphicsView::setCompositeMode(TomoGui_Constants::CompositeType mode)
{
  m_ImageDisplayType = TomoGui_Constants::CompositedImage;
  switch(mode)
  {
    case TomoGui_Constants::Exclusion: m_composition_mode = QPainter::CompositionMode_Exclusion; break;
    case TomoGui_Constants::Difference: m_composition_mode = QPainter::CompositionMode_Difference; break;
    case TomoGui_Constants::Alpha_Blend:
      m_composition_mode = QPainter::CompositionMode_SourceOver;
      break;
#if 0
    case 2: m_composition_mode = QPainter::CompositionMode_Plus; break;
    case 3: m_composition_mode = QPainter::CompositionMode_Multiply; break;
    case 4: m_composition_mode = QPainter::CompositionMode_Screen; break;
    case 5: m_composition_mode = QPainter::CompositionMode_Darken; break;
    case 6: m_composition_mode = QPainter::CompositionMode_Lighten; break;
    case 7: m_composition_mode = QPainter::CompositionMode_ColorDodge; break;
    case 8: m_composition_mode = QPainter::CompositionMode_ColorBurn; break;
    case 9: m_composition_mode = QPainter::CompositionMode_HardLight; break;
    case 10: m_composition_mode = QPainter::CompositionMode_SoftLight; break;

    case 12: m_composition_mode = QPainter::CompositionMode_Destination; break;
    case 13: m_composition_mode = QPainter::CompositionMode_Source; break;
    case 14: m_composition_mode = QPainter::CompositionMode_DestinationOver; break;
    case 15: m_composition_mode = QPainter::CompositionMode_SourceIn; break;
    case 16: m_composition_mode = QPainter::CompositionMode_DestinationIn; break;
    case 17: m_composition_mode = QPainter::CompositionMode_DestinationOut; break;
    case 18: m_composition_mode = QPainter::CompositionMode_SourceAtop; break;
    case 19: m_composition_mode = QPainter::CompositionMode_DestinationAtop; break;
    case 20: m_composition_mode = QPainter::CompositionMode_Overlay; break;
    case 21: m_composition_mode = QPainter::CompositionMode_Clear; break;
#endif
  default:
    m_composition_mode = QPainter::CompositionMode_Exclusion; break;
  }

  this->setImageDisplayType(m_ImageDisplayType);
}
